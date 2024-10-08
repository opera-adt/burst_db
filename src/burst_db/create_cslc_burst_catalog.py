import ast
import csv
import logging
import re
import tempfile
from dataclasses import astuple
from datetime import date, datetime
from pathlib import Path

import click
import duckdb
import pandas as pd
from opera_utils import missing_data
from tqdm.contrib.concurrent import thread_map

logger = logging.getLogger("burst_db")


def create_burst_catalog(input_csv: Path, opera_db: Path, output_file: Path):
    """Create a burst catalog from the input CSV and OPERA database.

    Parameters
    ----------
    input_csv : Path
        Path to the input CSV file containing burst data.
    opera_db : Path
        Path to the OPERA database containing the frames table.
    output_file : Path
        Path to the output file (can be .csv, .parquet, or .duckdb).

    """
    logging.info("Creating table from CMR raw bursts...")
    with duckdb.connect(str(output_file)) as conn:
        # Note: we need to add the `FIRST` because CMR has duplicates.
        # the csv has an odd header with a #
        # the "revision time" is sort of like "PGE Processing time"
        # The "Temporal Time" is the acquisition sensing time in the S1 filename
        conn.sql(
            f"""
    CREATE TABLE bursts AS
    SELECT
        LOWER(REPLACE(substring("# Granule ID", 18, 15), '-', '_')) AS burst_id_jpl,
        "Temporal Time"::TIMESTAMP AS sensing_time,
        MAX("Revision Time"::TIMESTAMP) AS max_revision_time,
        FIRST("# Granule ID") AS granule_id,
        granule_id[72:73] AS pol,
        FIRST("Revision-Temporal Delta Hours") AS delta_hours,
        FIRST("revision-id") AS revision_id
    FROM read_csv_auto('{input_csv}', header=True, sample_size=-1)
    GROUP BY
        1, 2
    ORDER BY
        1, 2
    """
        )

        conn.sql(f"ATTACH '{opera_db}' AS opera")

        logging.info(
            "Joining to OPERA Frame/Burst database Creating table from CMR raw"
            " bursts..."
        )
        conn.sql(
            """
    CREATE TABLE bursts_with_frame_ids AS
    SELECT
        b.burst_id_jpl,
        b.sensing_time,
        f.fid AS frame_id,
        f.is_north_america
    FROM bursts b
    JOIN opera.burst_id_map bm ON b.burst_id_jpl = bm.burst_id_jpl
    JOIN opera.frames_bursts fb ON fb.burst_ogc_fid = bm.OGC_FID
    JOIN opera.frames f ON fb.frame_fid = f.fid
    WHERE pol = 'VV'
    """
        )


def fetch_bursts(db_file: Path | str):
    """Fetch bursts from the created database.

    Parameters
    ----------
    db_file : Path
        Path to the database file containing the bursts_with_frame_ids table.

    Returns
    -------
    duckdb.DuckDBDataFrame
        DataFrame containing the fetched bursts.

    """
    query = """
    SELECT
      frame_id,
      burst_id_jpl,
      sensing_time
    FROM
      bursts_with_frame_ids
    WHERE
      is_north_america
    ORDER BY
      frame_id,
      burst_id_jpl,
      sensing_time
    """
    with duckdb.connect(str(db_file)) as conn:
        df_out = conn.sql(query).df()

    df_out["sensing_date"] = df_out["sensing_time"].dt.date
    return df_out


columns_full = [
    "frame_id",
    "option_num",
    "total_num_bursts",
    "burst_id_list",
    "date_list",
    "num_candidate_bursts",
    "inputs",
]


def _date_tup_string_to_list(input_string):
    date_pattern = r"datetime\.date\((\d{4}), (\d{1,2}), (\d{1,2})\)"
    matches = re.findall(date_pattern, input_string)

    # Convert the extracted date components to datetime.date objects
    return [date(int(year), int(month), int(day)) for year, month, day in matches]


def make_consistent_burst_json(
    db_file: Path | str, date_range: str | None = None
) -> Path:
    """Create the JSON database of consistent burst IDs over time.

    Parameters
    ----------
    db_file : Path | str
        Path to the duckdb database of CSLC archived bursts in CMR.
    date_range : str | None, optional
        String describing date range that `db_file` covers.
        Used to name the output JSON file
        By default None

    Returns
    -------
    Path
        Path to output JSON file.

    """
    df_bursts_all = fetch_bursts(db_file=db_file)
    frame_ids = df_bursts_all.frame_id.unique()
    frame_id_date_to_sensing_time = (
        # Take the sensing time of the earlier burst id
        df_bursts_all.groupby(["frame_id", "sensing_date"])["sensing_time"]
        # only use the earliest burst time (the beginning sensing of this frame
        .min()
        # truncate away the microseconds
        .dt.floor("s").to_dict()
    )
    with tempfile.TemporaryDirectory() as tmpdir:
        out_dir = Path(tmpdir)
        out_dir.mkdir(exist_ok=True, parents=True)

        def to_rows_full(frame_id: int, options):
            return [(frame_id, i, *astuple(b)) for i, b in enumerate(options)]

        def get_options(frame_id: int):
            tups = (
                df_bursts_all[df_bursts_all.frame_id == frame_id][
                    ["burst_id_jpl", "sensing_date"]
                ]
                .to_numpy()
                .tolist()
            )
            return missing_data.get_missing_data_options(burst_id_date_tuples=tups)

        def run_search(frame_id: int):
            options = get_options(frame_id)
            rows_full = list(to_rows_full(frame_id, options))
            cur_name = out_dir / f"full_frame_id_{frame_id:05d}.csv"
            with open(cur_name, "w", newline="") as f:
                writer = csv.writer(
                    f, delimiter=",", quotechar='"', quoting=csv.QUOTE_MINIMAL
                )
                writer.writerow(columns_full)
                for opt in rows_full:
                    writer.writerow(opt)

        _ = thread_map(
            run_search,
            frame_ids,
            max_workers=4,
            desc="Finding consistent bursts per frame",
        )
        query = f"""COPY (
            SELECT
                *
            FROM
                '{out_dir}/full_frame*.csv'
            ORDER BY
                2,
                1
            ) TO
            '{out_dir}/missing_data_per_frame_full.csv'
            (HEADER TRUE);"""
        with duckdb.connect(str(db_file)) as conn:
            conn.sql(query)

        df_selected_as_lists = pd.read_csv(f"{out_dir}/missing_data_per_frame_full.csv")

    # now work with the processed CSV
    df_selected_as_lists.columns = columns_full
    df_selected_as_lists = df_selected_as_lists.loc[
        df_selected_as_lists.option_num == 0
    ]
    df_selected_as_lists["burst_id_list"] = df_selected_as_lists["burst_id_list"].apply(
        lambda x: list(ast.literal_eval(x))
    )
    df_selected_as_lists["date_list"] = df_selected_as_lists["date_list"].apply(
        _date_tup_string_to_list
    )

    df_selected_as_lists["sensing_time_list"] = df_selected_as_lists.apply(
        lambda row: [
            frame_id_date_to_sensing_time.get((row["frame_id"], date), None).strftime(
                "%Y-%m-%dT%H:%M:%S"
            )
            for date in row["date_list"]
        ],
        axis=1,
    )
    today_str = datetime.today().strftime("%Y-%m-%d")
    date_str = today_str if date_range is None else date_range
    base_name = f"opera-disp-s1-consistent-burst-ids-{today_str}"
    if date_range is not None:
        base_name += f"-{date_str}"

    output_file = Path(base_name).with_suffix(".json")
    df_output = df_selected_as_lists.set_index("frame_id")[
        ["burst_id_list", "sensing_time_list"]
    ]
    df_output.to_json(output_file, orient="index", date_format="iso", indent=2)

    return output_file


@click.command()
@click.argument("input_csv", type=click.Path(exists=True, path_type=Path))
@click.argument("opera_db", type=click.Path(exists=True, path_type=Path))
@click.argument("output_file", type=click.Path(path_type=Path))
def make_burst_catalog(input_csv: Path, opera_db: Path, output_file: Path):
    """Create a burst catalog and consistent burst JSON.

    INPUT_CSV: Path to the input CSV file containing CMR-queried burst data.
    OPERA_DB: Path to the OPERA DISP Frame geopackage database.
    OUTPUT_FILE: Path to the output .duckdb file for the full deduped bursts.
    """
    if not output_file.exists():
        create_burst_catalog(input_csv, opera_db, output_file)
        logger.info(f"Burst catalog created: {output_file}")
    else:
        logger.info(f"Using existing burst catalog {output_file}")
    try:
        # Name should be like:
        # "cmr_survey.csv.raw.2016-07-01_to_2024-09-04.csv"
        # So grab the last part
        date_range = input_csv.stem.split(".")[-1]
    except Exception:
        logger.warning(f"Failed to parse date range from {input_csv}", exc_info=True)
        date_range = None

    make_consistent_burst_json(output_file, date_range=date_range)
