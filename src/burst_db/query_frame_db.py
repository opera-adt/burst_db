# /usr/bin/env python
import json
import sqlite3
from pathlib import Path
from typing import Optional, Tuple

import click
import pandas as pd

from .utils import build_wkt_from_bbox

DEFAULT_DB = Path("~/dev/opera-s1-disp.gpkg").expanduser()

Bbox = tuple[float, float, float, float]
MaybeBbox = Optional[Bbox]


def query_database(frame_id: int, db_path: Path) -> dict:
    """
    Query the geopackage database with the provided frame ID.

    Parameters
    ----------
    frame_id : int
        The frame ID to query.
    db_path : str
        Path to the geopackage database.

    Returns
    -------
    dict
        Result of the query.

    """
    query = """
            SELECT
                f.fid AS frame_id,
                f.epsg,
                f.is_land,
                f.is_north_america,
                MIN(xmin) AS xmin,
                MIN(ymin) AS ymin,
                MAX(xmax) AS xmax,
                MAX(ymax) AS ymax,
                GROUP_CONCAT(burst_id_jpl) AS burst_ids
            FROM frames f
            JOIN frames_bursts fb ON fb.frame_fid = f.fid
            JOIN burst_id_map b ON fb.burst_ogc_fid = b.ogc_fid
            WHERE f.fid = ?
            GROUP BY 1;
    """
    with sqlite3.connect(db_path) as con:
        df_frame_to_burst = pd.read_sql_query(query, con, params=(frame_id,))

    df_frame_to_burst.burst_ids = df_frame_to_burst.burst_ids.str.split(",")
    df_frame_to_burst.is_land = df_frame_to_burst.is_land.astype(bool)
    df_frame_to_burst.is_north_america = df_frame_to_burst.is_north_america.astype(bool)
    out_dict = df_frame_to_burst.set_index("frame_id").to_dict(orient="index")
    return out_dict[frame_id]


@click.command(context_settings={"show_default": True})
@click.option(
    "--db-path",
    type=click.Path(exists=True),
    default=DEFAULT_DB,
    show_default=True,
    help="Path to the geopackage database.",
)
@click.option(
    "--bbox",
    type=float,
    nargs=4,
    help="(left, bottom, right, top) Lon/lat bounding box of area of interest.",
)
@click.option("--wkt", help="Well-Known Text (WKT) representation of geometry.")
def intersect(
    db_path: Path, bbox: Optional[Tuple[float, float, float, float]], wkt: Optional[str]
):
    """Query for frames intersecting a given bounding box or WKT geometry."""
    if bbox is None and wkt is None:
        raise click.BadParameter(
            "Please provide either --bbox or --wkt option, not both or neither."
        )
    wkt_str = build_wkt_from_bbox(*bbox) if bbox is not None else wkt
    # Doing this query is slow, doesn't use the spatial index
    # query = """
    # SELECT *
    # FROM burst_id_map
    # WHERE Intersects(PolygonFromText(?, 4326), GeomFromGPB(geom))
    # """
    # For use of RTree:
    # https://erouault.blogspot.com/2017/03/dealing-with-huge-vector-geopackage.html
    # Also, need the `GeomFromGPB` to convert geopackage blob to spatialite formate
    # http://www.gaia-gis.it/gaia-sins/spatialite-sql-4.2.0.html
    query = """
WITH given_geom AS (
    SELECT PolygonFromText(?, 4326) as g
),
BBox AS (
    SELECT
        MbrMinX(g) AS minx,
        MbrMaxX(g) AS maxx,
        MbrMinY(g) AS miny,
        MbrMaxY(g) AS maxy
    FROM given_geom
)

SELECT fid as frame_id, epsg, relative_orbit_number, orbit_pass,
       is_land, is_north_america, ASText(GeomFromGPB(geom)) AS wkt
FROM frames
WHERE fid IN (
    SELECT id
    FROM rtree_frames_geom
    JOIN BBox ON
        rtree_frames_geom.minx <= BBox.maxx AND
        rtree_frames_geom.maxx >= BBox.minx AND
        rtree_frames_geom.miny <= BBox.maxy AND
        rtree_frames_geom.maxy >= BBox.miny
)
AND Intersects((SELECT g FROM given_geom), GeomFromGPB(geom));
"""

    with sqlite3.connect(db_path) as con:
        con.enable_load_extension(True)
        con.load_extension("mod_spatialite")
        df_intersecting_frames = pd.read_sql_query(query, con, params=[wkt_str])

    out_dict = df_intersecting_frames.to_dict(orient="records")
    click.echo(json.dumps(out_dict, indent=4))


@click.command(context_settings={"show_default": True})
@click.argument("frame_id", type=int)
@click.option(
    "--db-path",
    type=click.Path(exists=True),
    default=DEFAULT_DB,
    show_default=True,
    help="Path to the geopackage database.",
)
def lookup(frame_id: int, db_path: Path):
    """Query the geopackage database for one frame ID."""
    # Assuming `query_database` is a function defined elsewhere that returns the result
    result = query_database(frame_id, db_path)
    click.echo(json.dumps(result, indent=4))
