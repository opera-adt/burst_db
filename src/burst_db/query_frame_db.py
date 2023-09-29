# /usr/bin/env python
import json
import sqlite3
from pathlib import Path

from typing import Optional
from shapely import geometry, box
import pandas as pd
import typer
from typing_extensions import Annotated

DEFAULT_DB = Path("~/dev/opera-s1-disp.gpkg").expanduser()


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


def build_wkt_from_bbox(xmin: float, ymin: float, xmax: float, ymax: float) -> str:
    """Convert bounding box coordinates to WKT POLYGON string."""
    return box(xmin, ymin, xmax, ymax).wkt


def intersect(
    db_path: Path = DEFAULT_DB,
    bbox: Optional[tuple[float, float, float, float]] = typer.Option(
        None,
        "--bbox",
        metavar="LEFT BOTTOM RIGHT TOP",
        help="Bounding box in format 'xmin,ymin,xmax,ymax'.",
    ),
    wkt: Optional[str] = typer.Option(
        None, "--wkt", help="Well-Known Text (WKT) representation of geometry."
    ),
):
    """Query for frames intersecting a given bounding box or WKT geometry."""

    if not (bool(bbox) ^ bool(wkt)):
        raise typer.BadParameter(
            "Please provide either --bbox or --wkt option, not both or neither."
        )

    if bbox:
        wkt_str = build_wkt_from_bbox(*bbox)
    else:
        wkt_str = wkt

    query = """
    SELECT OGC_FID, burst_id_jpl, epsg, relative_orbit_number, orbit_pass, time_from_anx_sec, ASText(geom) AS wkt
    FROM burst_id_map
    WHERE Intersects(PolygonFromText(?, 4326), geom)
    """
    typer.echo((query, wkt_str))

    with sqlite3.connect(db_path) as con:
        con.enable_load_extension(True)
        con.load_extension("mod_spatialite")
        df_intersecting_frames = pd.read_sql_query(query, con, params=[wkt_str])

    print(df_intersecting_frames)


def lookup(
    frame_id: int,
    db_path: Annotated[
        Path, typer.Option(help="Path to the geopackage database.")
    ] = DEFAULT_DB,
    # db_path: Path = DEFAULT_DB,
):
    """Query the geopackage database and return the result as JSON based on the provided frame ID."""
    result = query_database(frame_id, db_path)
    typer.echo(json.dumps(result, indent=4))
