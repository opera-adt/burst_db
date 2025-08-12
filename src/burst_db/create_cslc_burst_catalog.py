import ast
import csv
import json
import logging
import re
import tempfile
from dataclasses import astuple
from datetime import date, datetime
from pathlib import Path
from typing import Optional

import click
import duckdb
import pandas as pd
from opera_utils import missing_data
from tqdm.contrib.concurrent import thread_map

logger = logging.getLogger("burst_db")

AUSTRALIAN_SAMPLE_FRAMES = (31577, 39355, 39360, 39362, 43655)
EDGE_FRAMES_TO_IGNORE = (
    28471,
    4854,
    12900,
    22674,
    40027,
    36261,
    42013,
)


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
        f.is_north_america as frame_is_north_america,
        bm.is_north_america as burst_is_north_america,
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
    query = f"""
    SELECT
      frame_id,
      burst_id_jpl,
      sensing_time
    FROM
      bursts_with_frame_ids
    WHERE
      (burst_is_north_america
      OR frame_id in {AUSTRALIAN_SAMPLE_FRAMES})
      AND frame_id not in {EDGE_FRAMES_TO_IGNORE}
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


def _date_is_excluded(
    frame_id: int | str, check_date: date, blackout_periods: dict[str, list[list[str]]]
) -> bool:
    if not blackout_periods.get(str(frame_id)):
        return False

    for start_str, end_str in blackout_periods[str(frame_id)]:
        start_date = datetime.fromisoformat(start_str).date()
        end_date = datetime.fromisoformat(end_str).date()
        if start_date <= check_date <= end_date:
            return True
    return False


def make_consistent_burst_json(
    db_file: Path | str,
    date_range: str | None = None,
    blackout_file: Path | str | None = None,
    input_files: Optional[dict[str, Path | str]] = None,
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
    blackout_file : Path | str, optional
        Filename of the per-frame blackout periods.
        If included, will skip over all sensing times a frame has within each period.
    input_files : dict[str, Path|str], optional
        For metadata tracking, the intput CMR survey/other files for trackign in
        the `metadata` key.

    Returns
    -------
    Path
        Path to output JSON file.

    """
    if input_files is None:
        input_files = {}
    df_bursts_all = fetch_bursts(db_file=db_file)
    if blackout_file:
        blackout_periods = json.loads(Path(blackout_file).read_text())["blackout_dates"]
    else:
        blackout_periods = {}

    frame_ids = df_bursts_all.frame_id.unique()
    frame_id_date_to_sensing_time = {
        (frame_id, sensing_date): sensing_time
        for (frame_id, sensing_date), sensing_time in (
            # Take the sensing time of the earlier burst id
            df_bursts_all.groupby(["frame_id", "sensing_date"])["sensing_time"]
            # only use the earliest burst time (the beginning sensing of this frame
            .min()
            # truncate away the microseconds
            .dt.floor("s").items()
        )
        # Filter away sensing times that fall in the blackout period
        if not _date_is_excluded(str(frame_id), sensing_date, blackout_periods)
    }
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

    def _to_dt_string(row):
        out = []
        for cur_date in row["date_list"]:
            d = frame_id_date_to_sensing_time.get((row["frame_id"], cur_date))
            if not d:
                continue
            out.append(d.strftime("%Y-%m-%dT%H:%M:%S"))
        return out

    df_selected_as_lists["sensing_time_list"] = df_selected_as_lists.apply(
        _to_dt_string, axis=1
    )

    today_str = datetime.today().strftime("%Y-%m-%d")
    basename = f"opera-disp-s1-consistent-burst-ids-{today_str}"
    if date_range is not None:
        basename = basename + f"-{date_range}"
    output_file = Path(basename).with_suffix(".json")

    df_output = df_selected_as_lists.set_index("frame_id")[
        ["burst_id_list", "sensing_time_list"]
    ]
    # Add "ignore_list" for time ranges to skip entirely
    # add to json, deliver
    out_data = df_output.to_dict(orient="index")
    total_out = {
        "metadata": {
            "generation_time": datetime.today(),
            "blackout_file": blackout_file,
            **input_files,
        },
        "data": out_data,
    }
    with open(output_file, "w") as f:
        f.write(json.dumps(total_out, indent=2, default=str))

    return output_file


@click.command()
@click.argument("input_csv", type=click.Path(exists=True, path_type=Path))
@click.argument("opera_db", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--full-db-path",
    type=click.Path(path_type=Path),
    help=(
        "Name of .duckdb file to store all burst catalog info. Default is"
        " `cslc-burst-database-{today}.duckdb`"
    ),
)
@click.option(
    "--blackout-file",
    type=click.Path(exists=True, path_type=Path),
    help="Path to per-frame blackout periods made by `create_blackout_dates_s1`",
)
def make_burst_catalog(
    input_csv: Path,
    opera_db: Path,
    full_db_path: Path | None,
    blackout_file: Path | None,
):
    """Create a burst catalog and consistent burst JSON.

    INPUT_CSV: Path to the input CSV file containing CMR-queried burst data.
    OPERA_DB: Path to the OPERA DISP Frame geopackage database.
    """
    if full_db_path is None:
        today = datetime.now().strftime("%Y-%m-%d")
        full_db_path = Path(f"cslc-burst-database-{today}.duckdb")

    if not full_db_path.exists():
        logger.info(f"Creating {full_db_path}")
        create_burst_catalog(input_csv, opera_db, full_db_path)
        logger.info(f"Burst catalog created: {full_db_path}")
    else:
        logger.info(f"Using existing burst catalog {full_db_path}")

    # Name should be like:
    # "cmr_survey.csv.raw.2016-07-01_to_2024-09-04.csv"
    # or cmr_survey_2016-07-01_to_2024-12-10.csv
    pattern = r"(\d{4}-\d{2}-\d{2}_to_\d{4}-\d{2}-\d{2})"
    if match := re.search(pattern, input_csv.stem):
        date_range = match.group(1)
    else:
        date_range = None

    make_consistent_burst_json(
        full_db_path,
        date_range=date_range,
        blackout_file=blackout_file,
        input_files={
            "input_cmr_csv": input_csv,
            "opera_db_file": opera_db,
            "generated_burst_db": full_db_path,
        },
    )
