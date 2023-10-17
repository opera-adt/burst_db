from __future__ import annotations

import csv
import sqlite3
import subprocess
from datetime import datetime
from pathlib import Path

import geopandas as gpd
import pandas as pd
from shapely import wkt
from tqdm.auto import tqdm

# find S1* -name "*csv.zip" -exec cat {} + > all_bursts_v1.csv
# time for f in `find . -name "*.csv.zip" | sort`; do
#     unzip -p $f | awk -F';' 'BEGIN {OFS=";"} \
#           {gsub(/\.SAFE$/, "", $4); print}' >> all_bursts.csv;
#     echo "Done with $f";
# done
# find missing_downloaded/csvs -name "*csv" -exec cat {} + > missing_downloaded2.csv
# cat all_bursts_v1.csv missing_downloaded.csv > combined.csv


# Define a function to process each chunk
def _process_chunk(chunk_df, output_gpkg: Path):
    chunk_df["geometry"] = chunk_df["geometry"].apply(wkt.loads)
    gdf = gpd.GeoDataFrame(chunk_df, geometry="geometry")

    # if file exists, append, else create new file
    mode = "a" if output_gpkg.exists() else "w"
    gdf.to_file(output_gpkg, layer="bursts", driver="GPKG", mode=mode)


def csv_to_gpkg(input_file: str | Path, chunk_size=100_000):
    # Loop through chunks
    output_gpkg = Path(input_file).with_suffix(".gpkg")
    for chunk in pd.read_csv(
        input_file,
        delimiter=";",
        chunksize=chunk_size,
        names=["burst_id_jpl", "sensing_time", "geometry", "granule"],
    ):
        _process_chunk(chunk, output_gpkg)


def _estimate_num_lines(csv_file):
    """Give a rough estimate of number of lines using the first CSV line."""
    with open(csv_file) as f:
        line1_length = len(next(f))
    size_in_bytes = Path(csv_file).stat().st_size
    return int(size_in_bytes / line1_length)


def split_csv_by_columns(input_file: str | Path, write_headers: bool = False):
    with open(input_file, "r") as infile:
        reader = csv.reader(infile, delimiter=";")

        # Open the output files
        with (
            open("burst_times.csv", "w", newline="") as times_file,
            open("burst_geometries.csv", "w", newline="") as geom_file,
            open("burst_granules.csv", "w", newline="") as granule_file,
        ):
            # Create CSV writers for the output files
            times_writer = csv.writer(times_file, delimiter=";")
            geom_writer = csv.writer(geom_file)
            granule_writer = csv.writer(granule_file)

            if write_headers:
                times_writer.writerow(["burst_id_jpl", "sensing_time"])
                geom_writer.writerow(["geometry"])
                granule_writer.writerow(["granule"])

            # Iterate over rows in the input file
            for row in tqdm(reader, total=_estimate_num_lines(input_file)):
                if len(row) != 4:  # Check if row has the expected number of columns
                    raise ValueError(f"Unexpected row: {row}")
                burst_id, sensing_time, geometry, granule = row

                # Write data to the output files
                times_writer.writerow([burst_id, sensing_time])
                geom_writer.writerow([geometry])
                granule_writer.writerow([granule])


def create_normalized_db(output_file: str | Path = "bursts.db"):
    # Load the data
    granules_df = pd.read_csv("burst_granules.csv", names=["name"])
    times_df = pd.read_csv(
        "burst_times.csv", sep=";", names=["burst_id_jpl", "sensing_time"]
    )

    # Create unique dataframes
    unique_bursts_df = (
        times_df[["burst_id_jpl"]]
        .drop_duplicates()
        .reset_index(drop=True)
        .reset_index()
        .rename(columns={"index": "pkey"})
    )
    unique_bursts_df.pkey += 1
    unique_granules_df = (
        granules_df.drop_duplicates()
        .reset_index(drop=True)
        .reset_index()
        .rename(columns={"index": "pkey"})
    )
    unique_granules_df.pkey += 1

    # Connect to the SQLite database
    conn = sqlite3.connect(output_file)
    cursor = conn.cursor()
    # Create granule table
    sql = """
CREATE TABLE granules (
    pkey INTEGER PRIMARY KEY AUTOINCREMENT,
    name TEXT NOT NULL UNIQUE
);
"""
    cursor.execute(sql)
    # Create burst table
    sql = """
CREATE TABLE burst_ids (
    pkey INTEGER PRIMARY KEY AUTOINCREMENT,
    burst_id_jpl TEXT NOT NULL UNIQUE
);
"""
    cursor.execute(sql)

    # Create mapping table
    sql = """
CREATE TABLE bursts (
    pkey INTEGER PRIMARY KEY AUTOINCREMENT,
    burst_pkey INTEGER,
    granule_pkey INTEGER,
    sensing_time TIMESTAMP NOT NULL,
    FOREIGN KEY (burst_pkey) REFERENCES burst_ids(pkey),
    FOREIGN KEY (granule_pkey) REFERENCES granules(pkey)
);
"""
    cursor.execute(sql)

    # Populate the tables
    unique_bursts_df.to_sql("burst_ids", conn, if_exists="append", index=False)
    unique_granules_df.to_sql("granules", conn, if_exists="append", index=False)

    # For the main 'bursts' table, you'll want to replace textual identifiers
    # with their corresponding foreign keys
    merged_df = pd.merge(times_df, unique_bursts_df, how="left", on="burst_id_jpl")
    merged_df.rename(columns={"pkey": "burst_pkey"}, inplace=True)

    merged_granule_df = pd.merge(granules_df, unique_granules_df, how="left", on="name")
    merged_df.loc[:, "granule_pkey"] = merged_granule_df["pkey"]

    merged_df = pd.merge(merged_df, unique_granules_df, how="left", on="name")
    merged_df.rename(columns={"pkey": "burst_pkey"}, inplace=True)

    final_df = merged_df.loc[:, ["burst_pkey", "granule_pkey", "sensing_time"]]
    final_df.to_sql("bursts", conn, if_exists="append", index=False)

    cursor.execute("CREATE INDEX idx_burst_id_jpl ON burst_ids (burst_id_jpl);")
    cursor.execute("CREATE INDEX idx_granule_name ON granules (name);")
    cursor.execute("CREATE INDEX idx_bursts_sensing_time on bursts(sensing_time);")


def denormalize_opera_frames(frames_db_filename: str, outfile: str):
    """Join the many-to-many relation in the OPERA frames DB into one table.

    Used to make simpler queries once we join to the historical bursts table.
    """
    query = """
    SELECT
        OGC_FID
        ,burst_id_jpl
        ,b.epsg as burst_epsg
        ,min(CAST (f.fid AS INTEGER)) as min_frame_id
        ,max(CAST (f.fid AS INTEGER)) as max_frame_id
        ,min(CAST (f.epsg AS INTEGER)) as min_frame_epsg
        ,max(CAST (f.epsg AS INTEGER)) as max_frame_epsg
        ,count(f.fid) as num_frames
        ,max(is_land) as is_land
        ,max(is_north_america) as is_north_america
    FROM frames f
    JOIN frames_bursts fb ON fb.frame_fid = f.fid
    JOIN burst_id_map b ON fb.burst_ogc_fid = b.ogc_fid
    GROUP BY 1, 2;"""
    _run_sqlite_output_csv(query, frames_db_filename, outfile)


def create_duckdb(
    filename: str = "all_bursts.duckdb", outfile_base: str = "all_bursts_with_is_land"
):
    """Combine the denormalized frame table and historical bursts.

    Exports results as parquet and CSV.
    Requires `duckdb`.
    """
    import duckdb

    con = duckdb.connect(filename)

    # Read in all the yearly burst CSVs
    query = """CREATE TABLE bursts AS
    (SELECT * FROM read_csv(
        'bursts_20*.csv',
        delim=';',
        header=false,
        columns={
            'burst_id_jpl': 'VARCHAR',
            'sensing_time': 'DATETIME',
            'geometry': 'VARCHAR',
            'granule': 'VARCHAR'
        }
    ));"""
    print(query)
    # con.sql(query)

    # Read in the denormalized opera-frame DB
    query = (
        "CREATE TABLE frame_burst_denormalized AS "
        "SELECT * FROM read_csv_auto('frame_bursts_denormalized.csv');"
    )
    print(query)
    con.sql(query)

    # Denormalize the historical bursts, append the metadata columns
    query = """CREATE TABLE bursts_denormalized AS
    SELECT
        b.burst_id_jpl
        ,sensing_time
        ,granule
        ,min_frame_id
        ,max_frame_id
        ,min_frame_epsg
        ,max_frame_epsg
        ,num_frames
        ,is_land
        ,is_north_america
    FROM bursts b
    JOIN frame_burst_denormalized fbd
    ON (b.burst_id_jpl = fbd.burst_id_jpl);"""
    print(query)
    con.sql(query)

    # Export the result to CSV and parquet
    query = f"COPY bursts_denormalized TO '{outfile_base}.csv' (FORMAT CSV);"
    print(query)
    con.sql(query)
    query = f"COPY bursts_denormalized TO '{outfile_base}.parquet' (FORMAT PARQUET);"
    print(query)
    con.sql(query)
    con.close()

    query = """CREATE TABLE bursts (
        burst_id_jpl VARCHAR,
        sensing_time TIMESTAMP,
        granule VARCHAR,
        min_frame_id  INTEGER,
        max_frame_id  INTEGER,
        min_frame_epsg  INTEGER,
        max_frame_epsg  INTEGER,
        num_frames INTEGER,
        is_land INTEGER,
        is_north_america INTEGER
    );
    """
    con2 = sqlite3.connect(f"{outfile_base}.sqlite3")
    con2.execute(query)

    query = f".import --skip 1 {outfile_base}.csv bursts"
    cmd = f'sqlite3 -init /dev/null -csv {outfile_base}.sqlite3 "{query}"'
    subprocess.run(cmd, shell=True)

    query = (
        "CREATE INDEX idx_burst_id_jpl_sensing_time on bursts(burst_id_jpl,"
        " sensing_time);"
    )
    print(query)
    con2.execute(query)


def _run_sqlite_output_csv(
    query: str,
    frames_db_filename: str,
    outfile: str,
):
    cmd = (
        f"sqlite3 -init /dev/null -header -csv {frames_db_filename} "
        f'"{query}"'
        f" > {outfile}"
    )
    print(cmd)

    subprocess.run(cmd, check=True, shell=True)


def make_first_seen_csv(
    frames_db_filename: str, outfile: str, minimum_date: datetime | None = None
):
    # QUERY: number of bursts per year
    # SELECT strftime('%Y', sensing_time) AS year_of_burst ...
    query = "SELECT burst_id_jpl, granule, min(sensing_time) as first_seen FROM bursts "
    if minimum_date is not None:
        query += f" WHERE sensing_time > '{minimum_date.strftime('%Y-%m-%d')}'"
    query += " GROUP BY 1"
    _run_sqlite_output_csv(query, frames_db_filename, outfile)


def get_bursts_over_na():
    from itertools import chain

    from burst_db import build_frame_db

    frame_to_burst = build_frame_db.read_zipped_json(
        "/home/staniewi/repos/burst_db/testframe/opera-s1-disp-frame-to-burst.json.zip"
    )
    opera_frames = [
        (k, v) for (k, v) in frame_to_burst["data"].items() if v["is_north_america"]
    ]
    opera_bursts = list(
        set(chain.from_iterable(v["burst_ids"] for k, v in opera_frames))
    )
    return opera_bursts, opera_frames


def normalize():
    all_csvs = list(Path(".").rglob("*.csv.zip"))

    def _read_bursts(f):
        return pd.read_csv(
            f,
            delimiter=";",
            names=["burst_id_jpl", "sensing_time", "geometry", "granule"],
        )

    df = pd.concat([_read_bursts(f) for f in all_csvs])
    df = pd.concat([_read_bursts(f) for f in all_csvs[:100]])

    df = pd.read_csv(
        "bursts_no_geom.csv",
        names=["burst_id_jpl", "sensing_time", "granule"],
        delimiter=";",
        parse_dates=["sensing_time"],
    )

    geoms = []
    with open("geoms.csv") as f:
        for line in f:
            geoms.append(wkt.loads(line))

    geoms = []
    ###
    idx = 0
    with open("bursts_geom_wkt.csv") as f:
        for line in f:
            geoms.append(wkt.loads(line))
            idx += 1
    print(f"Loaded {idx} lines")
    gdf = gpd.GeoDataFrame(df, geometry=geoms)

    gdf.set_crs(epsg=4326, inplace=True)
    gdf.to_file("all_bursts.gpkg", layer="bursts")

    gdf[["burst_id_jpl", "sensing_time", "geometry"]].to_file(
        "all_bursts_separated.gpkg", layer="bursts"
    )
    gdf[["granule"]].str.replace(".SAFE", "").to_file(
        "all_bursts_separated.gpkg", mode="a", layer="granules"
    )
    gdf.geometry.to_file("burst_geometries.gpkg", crs="EPSG:4326")
