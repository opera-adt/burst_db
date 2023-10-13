import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Optional

import typer
from disp_s1.utils import FRAME_TO_BURST_JSON_FILE, read_zipped_json

DB_PATH = Path("/home/staniewi/dev/coverage-map-s1/global/testing_all_bursts.gpkg")

from functools import partial

def fetch_granules(
    frame_ids: list[str],
    min_datetime: Optional[datetime] = None,
    max_datetime: Optional[datetime] = None,
    output_file: Optional[typer.FileTextWrite] = None,
    db_path: Path = DB_PATH,
    debug: bool = False,
) -> list[str]:
    # Step 1: Parse zipped JSON to get the burst_ids for the given frame_ids
    data = read_zipped_json(FRAME_TO_BURST_JSON_FILE)["data"]
    burst_ids = []
    for frame_id in frame_ids:
        burst_ids.extend(data.get(f"{int(frame_id)}", {}).get("burst_ids", []))
    # import ipdb
    # ipdb.set_trace()

    # Step 2: De-duplicate the burst_ids
    unique_burst_ids = list(set(burst_ids))

    # Step 3: Connect to the SQLite database and get the relevant granules
    with sqlite3.connect(db_path) as conn:
        if debug:
            # This adds a callback to print the executed statements to stderr
            conn.set_trace_callback(partial(typer.echo, err=True))

        cursor = conn.cursor()

        query = f"""
        SELECT DISTINCT granule
        FROM bursts
        WHERE bursts.burst_id_jpl IN ({', '.join(['?']*len(unique_burst_ids))})
        """
        args = unique_burst_ids
        if min_datetime:
            query += "\nAND sensing_time >= ?"
            args += [min_datetime.date()]

        if max_datetime:
            query += "\nAND sensing_time <= ?"
            args += [max_datetime.date()]

        results = cursor.execute(query, args).fetchall()

        # import ipdb 
        # ipdb.set_trace()
        for row in results:
            granule = row[0].replace(".SAFE", "")
            typer.echo(granule, file=output_file)
        