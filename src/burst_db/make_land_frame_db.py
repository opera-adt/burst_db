#!/usr/bin/env python
import argparse
import os
import shutil
import sqlite3
import subprocess
import tempfile
import time
import zipfile
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.affinity import translate
from shapely import STRtree

from burst_db import frames

ESA_DB_URL = "https://sar-mpc.eu/files/S1_burstid_20220530.zip"

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
    return output_path


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


def make_land_df(
    buffer_deg=0.2, outname="usgs_land_02deg_buffered.geojson", driver="GeoJSON"
):
    if outname and Path(outname).exists():
        return gpd.read_file(outname)
    df_land_cont = gpd.read_file("data/GSHHS_shp/h/GSHHS_h_L1.shp")
    df_antarctica = gpd.read_file("data/GSHHS_shp/h/GSHHS_h_L6.shp")
    df_land = pd.concat([df_land_cont, df_antarctica], axis=0)[["geometry"]].copy()
    df_land.geometry = df_land.geometry.buffer(buffer_deg)
    df_land = df_land.dissolve()
    if outname:
        df_land.to_file(outname, driver=driver)
    return df_land


def get_land_indicator(df_burst_triplet, land_geom):
    tree = STRtree(df_burst_triplet.geometry)
    idxs_land = tree.query(land_geom, predicate="intersects")
    if idxs_land.ndim == 2:
        idxs_land = idxs_land[1]
    is_in_land = df_burst_triplet.index.isin(idxs_land)
    return is_in_land


def create_frame_to_burst_mapping(is_in_land, min_frame, target_frame, max_frame):
    """Create the JOIN table between frames_number and burst_id."""
    indicator, consecutive_land_frames, land_slices = frames._buffer_small_frames(
        is_in_land, min_frame=min_frame
    )
    print("Number of occurrences with smallest consecutive land frames:")
    print(sorted(consecutive_land_frames.items())[:20])

    frame_slices = []
    for start_idx, end_idx in land_slices:
        cur_slices = frames.solve_dp(
            end_idx - start_idx,
            target=target_frame,
            min_frame=min_frame,
            max_frame=max_frame,
        )
        # bump up so they refer to rows, instead of being from 0
        cur_slices = [(s + start_idx, e + start_idx) for (s, e) in cur_slices]
        frame_slices.extend(cur_slices)

    # Create the frame IDs mapping to burst_id
    # (frame_id, OGC_FID)
    frame_ogc_fid_tuples = []
    for frame_id, (start_idx, end_idx) in enumerate(frame_slices, start=1):
        for burst_id in range(start_idx + 1, end_idx + 1):
            for ogc_fid in range(1 + 3 * (burst_id - 1), 4 + 3 * (burst_id - 1)):
                frame_ogc_fid_tuples.append((frame_id, ogc_fid))

    df_frame_to_burst_id = pd.DataFrame(
        frame_ogc_fid_tuples, columns=["frame_fid", "burst_ogc_fid"]
    )
    return df_frame_to_burst_id


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
            "CREATE INDEX IF NOT EXISTS idx_burst_id_map_OGC_FID ON burst_id_map (OGC_FID)"
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
        con.execute("UPDATE gpkg_geometry_columns SET geometry_type_name = 'MULTIPOLYGON';")


def get_epsg_codes(df):
    """Get the EPSG codes for all non-antimeridian polygons in a GeoDataFrame.

    Uses the UTM library to account for the oddities of the Zones near Norway [1]_.

    References
    ----------
    .. _[1]: http://www.jaworski.ca/utmzones.htm
    """
    import utm  # https://github.com/Turbo87/utm

    epsgs = np.zeros(len(df), dtype=int)

    # do the antimeridian frames first
    am_idxs = df.geometry.map(lambda geo: (geo.geom_type != "Polygon")).values
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


def add_gpkg_spatial_ref_sys(outfile, epsgs):
    unique_epsgs = np.unique(epsgs)
    with sqlite3.connect(outfile) as con:
        _setup_spatialite_con(con)
        sql = "SELECT gpkgInsertEpsgSRID({epsg});"
        for epsg in unique_epsgs:
            try:
                con.execute(sql.format(epsg=epsg))
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
    WITH bboxes(b, OGC_FID) AS (
        SELECT ST_Envelope(ST_Transform(geom, epsg)), OGC_FID
        FROM burst_id_map
        WHERE epsg IS NOT NULL 
        AND epsg != 0
    )
    UPDATE burst_id_map SET (xmin, ymin, xmax, ymax) = (
    SELECT
        FLOOR((ST_MinX(b) - {margin}) / {snap:.1f}) * {snap:.1f},
        FLOOR((ST_MinY(b) - {margin}) / {snap:.1f}) * {snap:.1f},
        CEIL((ST_MaxX(b) + {margin}) / {snap:.1f}) * {snap:.1f},
        CEIL((ST_MaxY(b) + {margin}) / {snap:.1f}) * {snap:.1f}
    FROM bboxes
    WHERE OGC_FID = bboxes.OGC_FID
    );
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
        con.execute("CREATE INDEX idx_burst_id_jpl on burst_id_map (burst_id_jpl);")


def get_parser():
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
        "--min-frame",
        type=int,
        default=5,
        help="Number of bursts per frame",
    )
    parser.add_argument(
        "--target-frame",
        type=int,
        default=10,
        help="Number of bursts per frame",
    )
    parser.add_argument(
        "--max-frame",
        type=int,
        default=12,
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
        default=0.2,
        help="If provided, a buffer (in degrees) to indicate that a frame is near land.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    t0 = time.time()

    args = get_parser()
    target_frame = args.target_frame
    min_frame = args.min_frame
    max_frame = args.max_frame
    if not args.outfile:
        basename = f"s1-frames-{target_frame}frames-{min_frame}min-{max_frame}max"
        outfile = f"{basename}.gpkg"

    esa_db_path = args.esa_db_path
    # Read ESA's Burst Data
    if not Path(esa_db_path).exists():
        print(f"Downloading {ESA_DB_URL} to {esa_db_path}...")
        esa_db_path = get_esa_burst_db()

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
    df_burst.to_file(outfile, driver="GPKG", layer="burst_id_map")

    print("Aggregating burst triplets (grouping IW1,2,3 geometries together)")
    df_burst_triplet = make_burst_triplets(df_burst)
    # Get the land polygon to intersect
    if args.land_buffer_deg is None:
        raise ValueError("Must provide a land buffer in degrees")

    print("Indicating which bursts are near land...")
    df_land = make_land_df(
        args.land_buffer_deg,
        outname=f"usgs_land_{args.land_buffer_deg}deg_buffered.geojson",
        driver="GeoJSON",
    )
    land_geom = df_land.geometry

    is_in_land = get_land_indicator(df_burst_triplet, land_geom)

    # Create frames and JOIN tables
    # Make the JOIN table first
    print("Defining frames - bursts JOIN table")
    df_frame_to_burst_id = create_frame_to_burst_mapping(
        is_in_land,
        min_frame=args.min_frame,
        target_frame=args.target_frame,
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
        df_frames = df_frames.loc[
            :, ["epsg", "geometry"]
        ]

    df_frames.to_file(outfile, driver="GPKG", layer="frames")

    update_burst_epsg(outfile)
    # Create the bounding box in UTM coordinates
    add_gpkg_spatial_ref_sys(outfile, epsgs)
    save_utm_bounding_boxes(outfile, margin=args.margin, snap=args.snap)

    # Make the minimal version of the DB
    ext = Path(outfile).suffix
    out_minimal = outfile.replace(ext, f"_bbox_only{ext}")
    print(f"Creating a epsg/bbox only version: {out_minimal}")
    make_minimal_db(outfile, out_minimal)

    tf = time.time()
    print(f"Total time: {tf - t0:.2f} seconds")
