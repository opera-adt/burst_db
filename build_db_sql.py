import argparse
import datetime
import os
import shutil
import sqlite3
import subprocess
import tempfile
import time
import zipfile
from multiprocessing import Pool, cpu_count

import pandas as pd
from shapely import ops, wkb

import utils as ut

ESA_DB_PATH = "burst_map_IW_000001_375887.sqlite3"
ESA_DB_URL = "https://sar-mpc.eu/files/S1_burstid_20220530.zip"


def convert_wkb_to_2d(
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

        t = time.time()
        df.rename(columns={"GEOMETRY": "geometry_wkb"}, inplace=True)
        geom_2d = _convert_wkb_to_2d_parallel(df["geometry_wkb"])
        df["geometry_wkb"] = geom_2d
        print(f"Converted {len(df)} rows in {time.time() - t:.2f} seconds")

    return df

def _convert_wkb_to_2d_parallel(wkb_series):
    workers = min(25, cpu_count())
    # return list(map(_convert_wkb_to_2d, wkb_series))
    print(f"Using {workers} workers to convert wkb to 2d")
    with Pool(workers) as p:
        # print(p.map(f, [1, 2, 3]))
        return list(p.map(ut._convert_wkb_to_2d, wkb_series))


def make_jpl_burst_db(df, db_path, table_name="burst_id_map"):
    with sqlite3.connect(db_path) as con:
        print("writing into database...")
        df.to_sql(table_name, con, if_exists="replace", index=False)
    sql = f"""
BEGIN;
SELECT InitSpatialMetaData();
-- Bug in spatialite: https://www.gaia-gis.it/fossil/libspatialite/tktview/8b6910dbbb2180026af54a5cc5aac107fb1d62ad

INSERT INTO spatial_ref_sys
(srid, auth_name, auth_srid, ref_sys_name, proj4text, srtext)
VALUES (32760, 'epsg', 32760,
'WGS 84 / UTM zone 60S',
'+proj=utm +zone=60 +south +datum=WGS84 +units=m +no_defs',
'PROJCS["WGS 84 / UTM zone 60S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32760"]]'
);

CREATE INDEX idx_OGC_FID on {table_name} (OGC_FID);

SELECT AddGeometryColumn('{table_name}', 'geometry', 4326, 'MULTIPOLYGON', 'XY');
UPDATE {table_name} SET geometry=GeomFromWKB(geometry_wkb, 4326);

ALTER TABLE {table_name} ADD epsg INTEGER;

-- Get the centroid of the geometry, use to find epsg code
WITH latlon AS
(
    SELECT Y(Centroid(geometry)) as lat,
           X(Centroid(geometry)) as lon,
           OGC_FID 
    FROM {table_name}
)
UPDATE {table_name}
SET epsg =
CASE
    WHEN lat >= 75   THEN 3413
    WHEN lat <= -60  THEN 3031
    WHEN lat > 0     THEN 32601 + CAST(ROUND((lon + 177) / 6, 0) AS INTEGER)
    ELSE                  32701 + CAST(ROUND((lon + 177) / 6, 0) AS INTEGER)
END
FROM latlon
WHERE {table_name}.OGC_FID = latlon.OGC_FID;

-- Columns for the bounding box limits in UTM
ALTER TABLE {table_name} ADD xmin FLOAT;
ALTER TABLE {table_name} ADD ymin FLOAT;
ALTER TABLE {table_name} ADD xmax FLOAT;
ALTER TABLE {table_name} ADD ymax FLOAT;

-- Get the bounding box and save the limits
WITH bboxes AS (
    SELECT envelope(transform(geometry, epsg)) as bb,
    OGC_FID
    FROM {table_name}
)
UPDATE {table_name} SET (xmin, ymin, xmax, ymax) = (
    Round((ST_MinX(bb) - 1000) / 50.0) * 50.0,
    Round((ST_MinY(bb) - 1000) / 50.0) * 50.0,
    Round((ST_MaxX(bb) + 1000) / 50.0) * 50.0,
    Round((ST_MaxY(bb) + 1000) / 50.0) * 50.0
)
FROM bboxes where bboxes.OGC_FID = {table_name}.OGC_FID;
;

-- Add the jpl burst id by concatenating the burst/track/subswath info
ALTER TABLE burst_id_map ADD COLUMN burst_id_jpl VARCHAR(16);
UPDATE {table_name} SET burst_id_jpl =
    't' || printf('%03d', relative_orbit_number) || '_'
               || printf('%06d', burst_id) || '_'
               || lower(subswath_name);

CREATE INDEX idx_burst_id_jpl on {table_name} (burst_id_jpl);

-- Drop unnecessary columns
ALTER TABLE {table_name} DROP COLUMN geometry_wkb;
-- These two are part of the burst_id_jpl
ALTER TABLE {table_name} DROP COLUMN burst_id;
ALTER TABLE {table_name} DROP COLUMN subswath_name;

COMMIT;
"""
    with sqlite3.connect(db_path) as con:
        con.enable_load_extension(True)
        con.load_extension("mod_spatialite")
        con.executescript(sql)
        # Clean up the database to reclaim space
        con.execute("VACUUM;")


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
        "--limit",
        type=int,
        help="For testing, limit the number of rows to process",
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    args = run_cli()
    if not os.path.exists(args.esa_db_path):
        print("ESA database not found, downloading...")
        t = time.time()
        get_esa_burst_db(args.esa_db_path)
        print(f"ESA database downloaded in {time.time() - t:.2f} seconds")
    t0 = time.time()
    df = convert_wkb_to_2d(
        esa_db_path=args.esa_db_path, limit=args.limit
    )
    print(f"Converted in {time.time() - t0:.2f} seconds")

    t1 = time.time()
    make_jpl_burst_db(df, args.output_path, table_name="burst_id_map")
    print(f"Created DB {args.output_path} in {time.time() - t1:.2f} seconds")

    print(f"Total script time: {time.time() - t0:.2f} seconds")
