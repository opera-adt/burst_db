import argparse
import datetime
import os
import shutil
import sqlite3
import subprocess
import tempfile
import time
import zipfile

import pandas as pd

import utils as ut

ESA_DB_PATH = "burst_map_IW_000001_375887.sqlite3"
ESA_DB_URL = "https://sar-mpc.eu/files/S1_burstid_20220530.zip"


def read_esa_db(
    esa_db_path=ESA_DB_PATH,
    limit=None,
):
    """Read in the ESA burst database and convert the wkb geometry to 2D."""
    print("Reading original ESA burst database")
    with sqlite3.connect(esa_db_path) as con:
        query = f"SELECT * FROM burst_id_map"
        if limit:
            query += f" LIMIT {limit}"
        print(query)

        t = time.time()
        df = pd.read_sql_query(query, con)
        print(f"Read {len(df)} rows in {time.time() - t:.2f} seconds")
        df.rename(columns={"GEOMETRY": "geometry_wkb"}, inplace=True)

    return df


def make_jpl_burst_db(
    df, db_path, table_name="burst_id_map", max_procs=None, snap=50.0, margin=4000.0
):
    print("Creating JPL burst ID")
    df["burst_id_jpl"] = ut.make_jpl_burst_id(df)
    max_procs = max_procs or 30

    t = time.time()
    proc_results = ut.process_geom_parallel(
        df["geometry_wkb"], max_procs=max_procs, snap=snap, margin=margin
    )
    df[["geometry_wkb", "epsg", "xmin", "ymin", "xmax", "ymax"]] = proc_results
    print(f"Processed {len(df)} geometries in {time.time() - t:.2f} seconds")

    with sqlite3.connect(db_path) as con:
        df.to_sql(table_name, con, if_exists="replace", index=False)

    print("Creating spatialite tables and geometry")
    sql = f"""
BEGIN;
SELECT InitSpatialMetaData();

CREATE INDEX idx_OGC_FID on {table_name} (OGC_FID);

SELECT AddGeometryColumn('{table_name}', 'geometry', 4326, 'MULTIPOLYGON', 'XY');
UPDATE {table_name} SET geometry=GeomFromWKB(geometry_wkb, 4326);
SELECT CreateSpatialIndex('burst_id_map', 'geometry');

CREATE INDEX idx_burst_id_jpl on {table_name} (burst_id_jpl);

-- Drop unnecessary blob column
ALTER TABLE {table_name} DROP COLUMN geometry_wkb;

COMMIT;
"""
    with sqlite3.connect(db_path) as con:
        con.enable_load_extension(True)
        con.load_extension("mod_spatialite")
        con.executescript(sql)
        # Clean up the database to reclaim space
        con.execute("VACUUM;")


def make_minimal_db(db_path, output_path):
    """Make a minimal database with only the burst_id_jpl, epsg, and bbox columns."""
    with sqlite3.connect(db_path) as con:
        df = pd.read_sql_query(
            "SELECT burst_id_jpl, epsg, xmin, ymin, xmax, ymax FROM burst_id_map", con
        )

    with sqlite3.connect(output_path) as con:
        df.to_sql("burst_id_map", con, if_exists="replace", index=False)
        con.execute("CREATE INDEX idx_burst_id_jpl on burst_id_map (burst_id_jpl);")


def get_esa_burst_db(output_path="esa_burst_map.sqlite3"):
    """Download the ESA burst database and convert to 2D."""
    # Download the ESA burst database

    print("Downloading ESA burst database")
    db_filename = "S1_burstid_20220530/IW/sqlite/burst_map_IW_000001_375887.sqlite3"
    cur_dir = os.getcwd()
    output_path = os.path.abspath(output_path)
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            os.chdir(tmpdir)
            subprocess.check_call(["wget", ESA_DB_URL])

            with zipfile.ZipFile(ESA_DB_URL.split("/")[-1], "r") as zip_ref:
                zip_ref.extract(db_filename)
                shutil.move(db_filename, output_path)
                shutil.rmtree(db_filename.split("/")[0])
        finally:
            os.chdir(cur_dir)


def run_cli():
    # Get the esa burst database path and output path
    now_str = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--esa-db-path",
        default="burst_map_IW_000001_375887.sqlite3",
        help=(
            "Path to the ESA sqlite burst database to convert, "
            f"downloaded from {ESA_DB_URL} . Will be downloaded if not exists."
        ),
    )
    parser.add_argument(
        "--output-path",
        type=str,
        default=f"burst_map_IW_000001_375887.OPERA-JPL.{now_str}.sqlite3",
        help="Path to the output database",
    )
    parser.add_argument(
        "--snap",
        type=float,
        default=50.0,
        help="Snap the bounding box to the nearest multiple of this value.",
    )
    parser.add_argument(
        "--margin",
        type=float,
        default=4000.0,
        help="Add this margin surrounding the bounding box of bursts.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        help="For testing, limit the number of rows to process",
    )
    parser.add_argument(
        "--max-procs",
        type=int,
        help="Max CPU count to use for processing geometries",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = run_cli()

    t0 = time.time()
    if not os.path.exists(args.esa_db_path):
        print("ESA database not found, downloading...")
        get_esa_burst_db(args.esa_db_path)
        print(f"ESA database downloaded in {time.time() - t0:.2f} seconds")
    df = read_esa_db(esa_db_path=args.esa_db_path, limit=args.limit)

    t1 = time.time()
    make_jpl_burst_db(
        df,
        args.output_path,
        table_name="burst_id_map",
        max_procs=args.max_procs,
        margin=args.margin,
        snap=args.snap,
    )
    print(f"Created DB {args.output_path} in {time.time() - t1:.2f} seconds")

    ext = os.path.splitext(args.output_path)[1]
    out_minimal = args.output_path.replace(ext, f"_bbox_only{ext}")
    print(f"Creating a epsg/bbox only version: {out_minimal}")
    make_minimal_db(args.output_path, out_minimal)

    print(f"Total script time: {time.time() - t0:.2f} seconds")
