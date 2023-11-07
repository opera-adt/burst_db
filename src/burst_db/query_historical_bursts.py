from __future__ import annotations

import csv
import sqlite3
from datetime import date, datetime
from functools import partial
from pathlib import Path
from typing import Any, Callable, Optional, Sequence

import click

from burst_db.utils import read_zipped_json  # Ensure this import is correct

FRAME_TO_BURST_JSON_FILE = Path(
    "/home/staniewi/dev/opera-s1-disp-frame-to-burst.json.zip"
)
DB_PATH = Path("/home/staniewi/dev/coverage-map-s1/global/testing_all_bursts.gpkg")


def _fetch_base(
    frame_ids: list[str],
    row_processor: Callable = lambda row: row,
    select_columns: Sequence[str] = ["*"],
    min_datetime: Optional[datetime] = None,
    max_datetime: Optional[datetime] = None,
    output_file: click.File = None,
    db_path: Path = DB_PATH,
    frame_to_burst_json_file: Path = FRAME_TO_BURST_JSON_FILE,
    headers: bool = True,
    debug: bool = False,
) -> list[Any]:
    # Step 1: Parse zipped JSON to get the burst_ids for the given frame_ids
    data = read_zipped_json(frame_to_burst_json_file)["data"]
    burst_ids = []
    for frame_id in frame_ids:
        burst_ids.extend(data.get(f"{int(frame_id)}", {}).get("burst_ids", []))

    # Step 2: De-duplicate the burst_ids
    unique_burst_ids = list(set(burst_ids))

    out = []

    # Step 3: Connect to the SQLite database and get the relevant granules
    with sqlite3.connect(db_path) as conn:
        if debug:
            # This adds a callback to print the executed statements to stderr
            conn.set_trace_callback(partial(click.echo, err=True))

        cursor = conn.cursor()

        query, args = _get_query(
            unique_burst_ids=unique_burst_ids,
            select_columns=select_columns,
            min_datetime=min_datetime,
            max_datetime=max_datetime,
        )
        results = cursor.execute(query, args).fetchall()
        writer = csv.writer(output_file)

        if headers:
            writer.writerow([c.replace("DISTINCT ", "") for c in select_columns])
        for row in results:
            out_row = row_processor(row)
            writer.writerow([out_row])
            out.append(out_row)
    return out


@click.command(context_settings={"show_default": True})
@click.argument("frame_ids", nargs=-1)
@click.option(
    "--min-datetime",
    type=click.DateTime(),
    help="Minimum datetime for filtering results.",
)
@click.option(
    "--max-datetime",
    type=click.DateTime(),
    help="Maximum datetime for filtering results.",
)
@click.option("--output-file", type=click.File("w"), help="Output file path.")
@click.option(
    "--db-path",
    type=click.Path(exists=True),
    default=DB_PATH,
    show_default=True,
    help="Database file path.",
)
@click.option(
    "--frame-to-burst-json-file",
    type=click.Path(exists=True),
    default=FRAME_TO_BURST_JSON_FILE,
    show_default=True,
    help="Frame to burst JSON file path.",
)
@click.option("--headers", is_flag=True, help="Include headers in the output.")
@click.option("--debug", is_flag=True, help="Enable debug mode.")
def fetch_granules(
    frame_ids,
    min_datetime,
    max_datetime,
    output_file,
    db_path,
    frame_to_burst_json_file,
    headers,
    debug,
):
    """Fetch granules for given frame IDs."""

    def row_processor(row):
        return row[0].replace(".SAFE", "")

    return _fetch_base(
        frame_ids=frame_ids,
        row_processor=row_processor,
        min_datetime=min_datetime,
        max_datetime=max_datetime,
        output_file=output_file,
        select_columns=["DISTINCT granule"],
        db_path=db_path,
        frame_to_burst_json_file=frame_to_burst_json_file,
        headers=headers,
        debug=debug,
    )


@click.command(context_settings={"show_default": True})
@click.argument("frame_ids", nargs=-1)
@click.option(
    "--min-datetime",
    type=click.DateTime(),
    help="Minimum datetime for filtering results.",
)
@click.option(
    "--max-datetime",
    type=click.DateTime(),
    help="Maximum datetime for filtering results.",
)
@click.option("--output-file", type=click.File("w"), help="Output file path.")
@click.option(
    "--db-path",
    type=click.Path(exists=True),
    default=DB_PATH,
    show_default=True,
    help="Database file path.",
)
@click.option(
    "--frame-to-burst-json-file",
    type=click.Path(exists=True),
    default=FRAME_TO_BURST_JSON_FILE,
    show_default=True,
    help="Frame to burst JSON file path.",
)
@click.option("--headers", is_flag=True, help="Include headers in the output.")
@click.option(
    "--with-granule", is_flag=True, help="Include granule information in the output."
)
@click.option("--debug", is_flag=True, help="Enable debug mode.")
def fetch_bursts(
    frame_ids,
    min_datetime,
    max_datetime,
    output_file,
    db_path,
    frame_to_burst_json_file,
    headers,
    with_granule,
    debug,
):
    """Get all (burst_id_jpl, sensing_time) for a list of frame ids."""
    select_columns = ["burst_id_jpl", "sensing_time"]
    if with_granule:
        select_columns.append("granule")
    return _fetch_base(
        frame_ids=frame_ids,
        # row_processor=row_processor,
        min_datetime=min_datetime,
        max_datetime=max_datetime,
        output_file=output_file,
        select_columns=select_columns,
        db_path=db_path,
        frame_to_burst_json_file=frame_to_burst_json_file,
        headers=headers,
        debug=debug,
    )


def _get_query(
    unique_burst_ids: Sequence[str],
    select_columns: Sequence[str],
    min_datetime: Optional[datetime] = None,
    max_datetime: Optional[datetime] = None,
) -> tuple[str, list]:
    """Build up the query to run on the full database."""
    # fid  geom  burst_id_jpl     sensing_time                granule
    query_base = f"""
    FROM bursts
    WHERE bursts.burst_id_jpl IN ({', '.join(['?']*len(unique_burst_ids))})"""

    query = f"SELECT {','.join(select_columns)}" + query_base
    args: list[str | date] = list(unique_burst_ids)
    if min_datetime:
        query += "\nAND sensing_time >= ?"
        args += [min_datetime.date()]

    if max_datetime:
        query += "\nAND sensing_time <= ?"
        args += [max_datetime.date()]
    return query, args


def get_gdf(
    wkt: str,
    filename: str = "./all_bursts_with_opera_db.duckdb",
    min_datetime: Optional[datetime] = None,
    max_datetime: Optional[datetime] = None,
    to_geodataframe: bool = False,
    debug: bool = False,
):
    import duckdb
    import shapely.wkt

    sql = f"""
    SELECT b.burst_id_jpl, sensing_time, ST_AsText(bm.geom) as geom_wkt
    FROM bursts b
    JOIN burst_id_map bm using (burst_id_jpl)
    WHERE ST_Intersects(bm.geom, ST_Extent(ST_GeomFromText('{wkt}')))
    """
    if min_datetime:
        # allow date or datetime
        d = min_datetime.date() if isinstance(min_datetime, datetime) else min_datetime
        sql += f"AND sensing_time >= '{d}'\n"
    if max_datetime:
        d = max_datetime.date() if isinstance(max_datetime, datetime) else max_datetime
        sql += f"AND sensing_time <= '{d}'\n"

    if debug:
        print(sql)
    with duckdb.connect(filename) as con:
        con.install_extension("spatial")
        con.load_extension("spatial")
        df = con.sql(sql).df()
    if to_geodataframe:
        import geopandas as gpd

        gdf = gpd.GeoDataFrame(
            df, crs="EPSG:4326", geometry=gpd.GeoSeries.from_wkt(df.geom_wkt)
        )
        return gdf.drop("geom_wkt", axis="columns")
    else:
        geom = df.geom_wkt.apply(shapely.wkt.loads)
        df["geometry"] = geom
        df.drop("geom_wkt", axis="columns", inplace=True)
        return df
