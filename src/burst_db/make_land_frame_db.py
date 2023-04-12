#!/usr/bin/env python
import argparse
import datetime
import sqlite3
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
import utm  # https://github.com/Turbo87/utm
from shapely import STRtree
from shapely.affinity import translate
from tqdm.auto import tqdm

from burst_db import __version__

from . import frames
from ._esa_burst_db import ESA_DB_URL, get_esa_burst_db
from ._land_usgs import get_land_df


def make_jpl_burst_id(df):
    burst_id_jpl = (
        "t"
        + df["relative_orbit_number"].astype(str).str.zfill(3)
        + "_"
        + df["burst_id"].astype(str).str.zfill(6)
        + "_"
        + df["subswath_name"].str.lower()
    )
    return burst_id_jpl


def _setup_spatialite_con(con):
    con.enable_load_extension(True)
    # Try the two versions, mac and linux with .so
    try:
        con.load_extension("mod_spatialite")
    except:
        con.load_extension("mod_spatialite.so")
    # Allow us to use spatialite functions in on GPKG files.
    # https://medium.com/@joelmalone/sqlite3-spatialite-and-geopackages-66a08485da6c
    con.execute("SELECT EnableGpkgAmphibiousMode();")


def make_burst_triplets(df_burst):
    def join_track_numbers(orbits: list) -> str:
        orbits = list(set(orbits))
        orbits_str = list(map(str, orbits))
        return ",".join(orbits_str)

    df_burst_triplet_temp = df_burst.dissolve(
        by="burst_id",
        aggfunc={
            "OGC_FID": ["min", "max"],
            "relative_orbit_number": join_track_numbers,
            "orbit_pass": "first",
        },
        as_index=False,
    )
    df_burst_triplet = df_burst_triplet_temp.reset_index(drop=True)
    df_burst_triplet.columns = [
        "burst_id",
        "geom",
        "OGC_FID_min",
        "OGC_FID_max",
        "relative_orbit_numbers",
        "look_direction",
    ]
    return df_burst_triplet


def get_land_indicator(df_burst_triplet: pd.DataFrame, land_geom):
    tree = STRtree(df_burst_triplet.geometry)
    idxs_land = tree.query(land_geom, predicate="intersects")
    if idxs_land.ndim == 2:
        idxs_land = idxs_land[1]
    is_in_land = df_burst_triplet.index.isin(idxs_land)
    return is_in_land


def make_frame_to_burst_table(outfile, df_frame_to_burst_id):
    with sqlite3.connect(outfile) as con:
        _setup_spatialite_con(con)

        df_frame_to_burst_id.to_sql("frames_bursts", con)
        con.execute(
            "CREATE INDEX IF NOT EXISTS idx_frames_bursts_burst_ogc_fid ON frames_bursts (burst_ogc_fid)"
        )
        con.execute(
            "CREATE INDEX IF NOT EXISTS idx_frames_bursts_frame_fid ON frames_bursts (frame_fid)"
        )
        con.execute(
            "CREATE INDEX IF NOT EXISTS idx_burst_id_map_burst_id ON burst_id_map (burst_id)"
        )


def make_frame_table(outfile):
    with sqlite3.connect(outfile) as con:
        _setup_spatialite_con(con)
        con.execute("CREATE TABLE frames " "(fid INTEGER PRIMARY KEY, epsg INTEGER)")
        # https://groups.google.com/g/spatialite-users/c/XcWvAk7vg0c
        # should add geom after the table is created
        # table_name , geometry_column_name , geometry_type , with_z , with_m , srs_id
        con.execute(
            "SELECT gpkgAddGeometryColumn('frames', 'geom', 'MULTIPOLYGON', 2, 0, 4326);"
        )

        # Aggregates burst geometries for each frame into one
        con.execute(
            """INSERT INTO frames(fid, geom)
            SELECT fb.frame_fid as fid,
                    ST_UnaryUnion(ST_Collect(geom)) as geom
            FROM burst_id_map b
            JOIN
                frames_bursts fb
                ON b.ogc_fid = fb.burst_ogc_fid
            GROUP BY 1;
        """
        )
        print("Creating indexes and spatial index...")
        con.execute("CREATE INDEX IF NOT EXISTS idx_frames_fid ON frames (fid)")
        con.execute("SELECT gpkgAddSpatialIndex('frames', 'geom') ;")
        # Extra thing so that QGIS recognizes "frames" better
        con.execute(
            "UPDATE gpkg_geometry_columns SET geometry_type_name = 'MULTIPOLYGON';"
        )
        # Set the relative_orbit_number as the most common value for each frame
        con.execute("ALTER TABLE frames ADD COLUMN relative_orbit_number INTEGER;")
        con.execute(
            """WITH frame_tracks AS (
                SELECT f.fid,
                    CAST(ROUND(AVG(b.relative_orbit_number)) AS INTEGER) AS relative_orbit_number
                FROM frames f
                JOIN frames_bursts fb ON f.fid = fb.frame_fid
                JOIN burst_id_map b ON fb.burst_ogc_fid = b.ogc_fid
                GROUP BY 1
            ) UPDATE frames SET relative_orbit_number = frame_tracks.relative_orbit_number
            FROM frame_tracks
            WHERE frames.fid = frame_tracks.fid;
            """
        )
        # Set the orbit_pass as the value from the first burst
        con.execute("ALTER TABLE frames ADD COLUMN orbit_pass TEXT;")
        con.execute(
            """WITH op AS
                (SELECT f.fid,
                        b.orbit_pass
                FROM frames f
                JOIN frames_bursts fb ON f.fid = fb.frame_fid
                JOIN burst_id_map b ON fb.burst_ogc_fid = b.ogc_fid),
            frame_orbits AS (
                SELECT fid,
                    FIRST_VALUE(orbit_pass) OVER (PARTITION BY fid) AS orbit_pass
                FROM op
                GROUP BY fid)
            UPDATE frames SET orbit_pass = frame_orbits.orbit_pass
            FROM frame_orbits
            WHERE frames.fid = frame_orbits.fid;
            """
        )


def get_epsg_codes(df):
    """Get the EPSG codes for all non-antimeridian polygons in a GeoDataFrame.

    Uses the UTM library to account for the oddities of the Zones near Norway [1]_.

    References
    ----------
    .. _[1]: http://www.jaworski.ca/utmzones.htm
    """

    def _is_on_antimeridian(geom):
        return geom.geom_type == "MultiPolygon" and len(geom.geoms) > 1

    epsgs = np.zeros(len(df), dtype=int)

    # do the antimeridian frames first
    am_idxs = df.geometry.map(_is_on_antimeridian).values
    epsgs[am_idxs] = df[am_idxs].geometry.map(antimeridian_epsg)

    # everything else
    # get the x, y (lon, lat) coords of all other rows
    other_coords = np.array(
        df[~am_idxs].geometry.map(lambda g: tuple(g.centroid.coords)[0]).tolist()
    )
    xs, ys = other_coords.T
    ys_full_size = np.ones(len(epsgs)) * np.nan
    ys_full_size[~am_idxs] = ys

    idxs = np.logical_and.reduce((~am_idxs, ys_full_size > 84))
    epsgs[idxs] = 3413

    idxs = np.logical_and.reduce((~am_idxs, ys_full_size < -60))
    epsgs[idxs] = 3031

    utm_idxs = np.logical_and(ys < 84, ys > -60)
    north_idxs = ys[utm_idxs] > 0

    # North hemisphere

    zones_north = [
        utm.from_latlon(y, x)[2]
        for (y, x) in zip(ys[utm_idxs][north_idxs], xs[utm_idxs][north_idxs])
    ]
    idxs = np.logical_and.reduce((~am_idxs, ys_full_size < 84, ys_full_size > 0))
    epsgs[idxs] = 32600 + np.array(zones_north)

    # South hemisphere

    zones_south = [
        utm.from_latlon(y, x)[2]
        for (y, x) in zip(ys[utm_idxs][~north_idxs], xs[utm_idxs][~north_idxs])
    ]
    idxs = np.logical_and.reduce((~am_idxs, ys_full_size > -60, ys_full_size < 0))
    epsgs[idxs] = 32700 + np.array(zones_south)
    return epsgs


def antimeridian_epsg(mp):
    """Calculate the EPSG of multipolygons along the antimeridian.

    Parameters
    ----------
    mp : shapely.geometry.MultiPolygon
        The multipolygon to calculate the EPSG for.

    Returns
    -------
    epsg : int
        The EPSG code for the multipolygon.

    Notes
    -----

    The EPSG code is calculated by taking the weighted average of the centroid of the
    polygons in the multipolygon. The centroid is weighted by the area of the polygon.
    The centroid is shifted by 360 degrees if it is in the western hemisphere.
    """
    y_c = mp.centroid.y
    # check north/south pole cases
    if y_c >= 84.0:
        return 3413
    elif y_c <= -60.0:
        return 3031

    # otherwise, do the weighted average of the shifted polygons to get the centroid
    A = 0
    x_weighted = 0
    # might have 2 or 3 polygons
    for g in mp.geoms:
        A += g.area
        if g.centroid.x < 0:
            g_shifted = translate(g, xoff=360)
            x_weighted += g_shifted.centroid.x * g.area
        else:
            x_weighted += g.centroid.x * g.area
    x_c = x_weighted / A

    base = 32600 if y_c > 0 else 32700
    # longitude 179 gets 32660 north of the equator
    # longitude -179 gets 32601
    zone_addition = 1 if x_c < 180 else 60
    return base + zone_addition


def update_burst_epsg(outfile):
    with sqlite3.connect(outfile) as con:
        _setup_spatialite_con(con)
        # add index
        con.execute(
            "CREATE INDEX IF NOT EXISTS idx_burst_id_map_burst_id_jpl ON burst_id_map (burst_id_jpl)"
        )
        # Set the EPSG on every burst from the frames
        print("Updating burst EPSGs to match frames...")
        sql = """WITH burst_epsgs AS (
                    SELECT b.OGC_FID,
                        f.epsg
                    FROM burst_id_map b
                    JOIN frames_bursts fb
                    ON b.OGC_FID = fb.burst_ogc_fid
                    JOIN frames f
                    ON fb.frame_fid = f.fid
                )
                UPDATE burst_id_map
                SET epsg = burst_epsgs.epsg
                FROM burst_epsgs
                WHERE burst_id_map.OGC_FID = burst_epsgs.OGC_FID;
        """
        con.execute(sql)


def add_gpkg_spatial_ref_sys(outfile):
    epsgs = [3031, 3413, 4326] + list(range(32601, 32661)) + list(range(32701, 32761))
    with sqlite3.connect(outfile) as con:
        _setup_spatialite_con(con)
        sql = "SELECT gpkgInsertEpsgSRID({epsg});"
        for epsg in tqdm(epsgs):
            try:
                con.execute(sql.format(epsg=epsg))
            except (sqlite3.OperationalError, sqlite3.IntegrityError):
                # exists
                pass
        # Fix the gpkg_spatial_ref_sys table for missing UTM zone 32760
        # https://www.gaia-gis.it/fossil/libspatialite/tktview/8b6910dbbb2180026af54a5cc5aac107fb1d62ad?plaintext
        sql = """INSERT INTO gpkg_spatial_ref_sys (
                    srs_name,
                    srs_id,
                    organization,
                    organization_coordsys_id,
                    definition
                    )
                VALUES (
                    'WGS 84 / UTM zone 60S',
                    32760,
                    'EPSG',
                    32760,
'PROJCS["WGS 84 / UTM zone 60S",GEOGCS["WGS 84",DATUM["WGS_1984",SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]],AUTHORITY["EPSG","6326"]],PRIMEM["Greenwich",0,AUTHORITY["EPSG","8901"]],UNIT["degree",0.0174532925199433,AUTHORITY["EPSG","9122"]],AUTHORITY["EPSG","4326"]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",177],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",10000000],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AXIS["Easting",EAST],AXIS["Northing",NORTH],AUTHORITY["EPSG","32760"]]'
                );
        """
        try:
            con.execute(sql)
        except (sqlite3.OperationalError, sqlite3.IntegrityError):
            # exists
            pass
        # More the entries from gpkg_spatial_ref_sys to spatial_ref_sys
        # so we can use the `ST_Transform` function
        con.execute("DROP TABLE IF EXISTS spatial_ref_sys;")
        sql = """
        CREATE TABLE spatial_ref_sys (
            srid       INTEGER NOT NULL PRIMARY KEY,
            auth_name  VARCHAR(256),
            auth_srid  INTEGER,
            srtext     VARCHAR(2048),
            proj4text  VARCHAR(2048)
        );"""
        con.execute(sql)
        sql = """
        INSERT INTO spatial_ref_sys
        SELECT
            srs_id AS srid,
            organization AS auth_name,
            organization_coordsys_id AS auth_srid,
            definition AS srtext,
            NULL
        FROM gpkg_spatial_ref_sys;
        """
        con.execute(sql)


def save_utm_bounding_boxes(outfile, margin=4000, snap=50.0):
    try:
        with sqlite3.connect(outfile) as con:
            for col in ["xmin", "ymin", "xmax", "ymax"]:
                con.execute(f"ALTER TABLE burst_id_map ADD COLUMN {col} INTEGER;")
    except sqlite3.OperationalError:
        # Already exists
        pass

    sql = f"""
WITH transformed(g, OGC_FID) AS
  (SELECT ST_Envelope(ST_Transform(geom, epsg)) g,
          OGC_FID
   FROM burst_id_map
   WHERE epsg != 0 )
UPDATE burst_id_map
SET (xmin,
     ymin,
     xmax,
     ymax) = (bboxes.xmin,
              bboxes.ymin,
              bboxes.xmax,
              bboxes.ymax)
FROM
  (SELECT OGC_FID,
          FLOOR((ST_MinX(g) - {margin}) / {snap:.1f}) * {snap:.1f} AS xmin,
          FLOOR((ST_MinY(g) - {margin}) / {snap:.1f}) * {snap:.1f} AS ymin,
          CEIL((ST_MaxX(g) + {margin}) / {snap:.1f}) * {snap:.1f} AS xmax,
          CEIL((ST_MaxY(g) + {margin}) / {snap:.1f}) * {snap:.1f} AS ymax
   FROM transformed) AS bboxes
WHERE burst_id_map.OGC_FID = bboxes.OGC_FID ;
    """
    with sqlite3.connect(outfile) as con:
        _setup_spatialite_con(con)
        con.execute(sql)


def make_minimal_db(db_path, output_path):
    """Make a minimal database with only the burst_id_jpl, epsg, and bbox columns."""
    with sqlite3.connect(db_path) as con:
        df = pd.read_sql_query(
            "SELECT burst_id_jpl, epsg, xmin, ymin, xmax, ymax FROM burst_id_map", con
        )
    # Make sure snapped coordinates as integers (~40% smaller than REAL)
    df["xmin"] = df["xmin"].astype(int)
    df["ymin"] = df["ymin"].astype(int)
    df["xmax"] = df["xmax"].astype(int)
    df["ymax"] = df["ymax"].astype(int)

    with sqlite3.connect(output_path) as con:
        df.to_sql("burst_id_map", con, if_exists="replace", index=False)
        # Skip making the index since we don't need super fast queries
        # con.execute("CREATE INDEX idx_burst_id_jpl on burst_id_map (burst_id_jpl);")


def create_metadata_table(db_path, args):
    """Make the metadata table with the arguments used to create the database."""
    df = pd.DataFrame(
        [
            {
                "margin": args.margin,
                "snap": args.snap,
                "min_frame": args.min_frame,
                "target_frame": args.target_frame,
                "max_frame": args.max_frame,
                "version": __version__,
                "land_buffer_deg": args.land_buffer_deg,
                "last_modified": datetime.datetime.now().isoformat(),
            }
        ]
    )
    with sqlite3.connect(db_path) as con:
        df.to_sql("metadata", con, if_exists="replace", index=False)


def get_cli_args():
    parser = argparse.ArgumentParser(
        description="Generate frames for Sentinel-1 data",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "--snap",
        type=float,
        default=30.0,
        help="Snap the bounding box to the nearest multiple of this value.",
    )
    parser.add_argument(
        "--margin",
        type=float,
        default=5000.0,
        help="Add this margin surrounding the bounding box of bursts.",
    )
    parser.add_argument(
        "--min-frame",
        type=int,
        default=5,
        help="Number of bursts per frame",
    )
    parser.add_argument(
        "--target-frame",
        type=int,
        default=9,
        help="Number of bursts per frame",
    )
    parser.add_argument(
        "--max-frame",
        type=int,
        default=10,
        help="Number of bursts per frame",
    )
    parser.add_argument(
        "-o",
        "--outfile",
        help="Output file name (default is "
        "'s1-frames-{target_frame}frames-{min_frame}min-{max_frame}max.gpkg'",
    )
    parser.add_argument(
        "--land-buffer-deg",
        type=float,
        default=0.3,
        help="A buffer (in degrees) to indicate that a frame is near land.",
    )

    return parser.parse_args()


def main():
    args = get_cli_args()

    t0 = time.time()
    target_frame = args.target_frame
    min_frame = args.min_frame
    max_frame = args.max_frame
    if not args.outfile:
        basename = f"s1-frames-{target_frame}frames-{min_frame}min-{max_frame}max"
        outfile = f"{basename}.gpkg"
    else:
        outfile = args.outfile

    esa_db_path = args.esa_db_path
    # Read ESA's Burst Data
    if not Path(esa_db_path).exists():
        print(f"Downloading {ESA_DB_URL} to {esa_db_path}...")
        get_esa_burst_db(esa_db_path)

    print("Loading burst data...")
    sql = "SELECT * FROM burst_id_map"
    with sqlite3.connect(esa_db_path) as con:
        df_burst = gpd.GeoDataFrame.from_postgis(
            sql, con, geom_col="GEOMETRY", crs="EPSG:4326"
        ).rename_geometry("geom")

    print("Forming string JPL id")
    jpl_ids = make_jpl_burst_id(df_burst)
    df_burst.loc[:, "burst_id_jpl"] = jpl_ids
    # placeholder to compute later
    df_burst.loc[:, "epsg"] = 0

    # Start the outfile with the ESA database contents
    print("Saving initial version")
    df_burst.set_index("OGC_FID").to_file(
        outfile, driver="GPKG", layer="burst_id_map", index=False
    )
    # Adjust the primary key so it still matches original OGC_FID
    with sqlite3.connect(outfile) as con:
        con.execute("ALTER TABLE burst_id_map RENAME COLUMN fid TO OGC_FID;")

    print("Aggregating burst triplets (grouping IW1,2,3 geometries together)")
    df_burst_triplet = make_burst_triplets(df_burst)
    # Get the land polygon to intersect
    if args.land_buffer_deg is None:
        raise ValueError("Must provide a land buffer in degrees")

    print("Indicating which bursts are near land...")
    df_land = get_land_df(args.land_buffer_deg)
    land_geom = df_land.geometry

    is_in_land = get_land_indicator(df_burst_triplet, land_geom)

    # Create frames and JOIN tables
    # Make the JOIN table first
    print("Defining frames - bursts JOIN table")
    df_frame_to_burst_id = frames.create_frame_to_burst_mapping(
        is_in_land,
        target_frame=args.target_frame,
        min_frame=args.min_frame,
        max_frame=args.max_frame,
    )
    make_frame_to_burst_table(outfile, df_frame_to_burst_id)

    # make the "frames" table
    print("Making frames table by aggregating burst geometries...")
    make_frame_table(outfile)
    df_frames = gpd.read_file(outfile, layer="frames")

    print("Computing EPSG codes for each frame...")
    epsgs = get_epsg_codes(df_frames)
    df_frames.loc[:, "epsg"] = epsgs

    print("Final number of frames:", len(df_frames))
    print("Saving frames...")
    # weird geopandas thing?
    if "geom" in df_frames.columns and df_frames._geometry_column_name == "geometry":
        df_frames = df_frames.loc[:, ["epsg", "geometry"]]

    df_frames.to_file(outfile, driver="GPKG", layer="frames")

    update_burst_epsg(outfile)

    # Create the bounding box in UTM coordinates
    add_gpkg_spatial_ref_sys(outfile)
    save_utm_bounding_boxes(outfile, margin=args.margin, snap=args.snap)

    # Make the minimal version of the DB
    ext = Path(outfile).suffix
    out_minimal = outfile.replace(ext, f"-bbox-only{ext}")
    print(f"Creating a epsg/bbox only version: {out_minimal}")
    make_minimal_db(outfile, out_minimal)

    # Add metadata to each
    create_metadata_table(outfile, args)
    create_metadata_table(out_minimal, args)

    print(f"Total time: {time.time() - t0:.2f} seconds")


if __name__ == "__main__":
    main()
