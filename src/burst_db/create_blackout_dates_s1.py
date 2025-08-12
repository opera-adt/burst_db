import json
import logging
from datetime import datetime, timedelta
from pathlib import Path

import click
import geopandas as gpd
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def _yearly_windows(
    start_ts: pd.Timestamp,
    end_ts: pd.Timestamp,
    years: range,
) -> list[list[str]]:
    """Generate blackout windows for every year in *years*.

    If the end month/day occurs *earlier* in the calendar than the start
    month/day (e.g. Nov-01 ➜ May-31), the end year is `year + 1`.
    """
    s_month, s_day = start_ts.month, start_ts.day
    e_month, e_day = end_ts.month, end_ts.day

    windows: list[list[str]] = []
    for yr in years:
        start = pd.Timestamp(year=yr, month=s_month, day=s_day)
        end_year = yr + (e_month < s_month or (e_month == s_month and e_day < s_day))
        end = pd.Timestamp(year=end_year, month=e_month, day=e_day) + pd.Timedelta(
            hours=23, minutes=59, seconds=59
        )
        windows.append([start.isoformat(), end.isoformat()])
    return windows


def snow_months_to_blackout_json(
    input_file: Path | str,
    max_default_duration: float = 180,
) -> dict:
    """Create a JSON of blackout periods for DISP-S1 from `snow_month_filter` outputs.

    The input dataframe has the follow fields for each row ():
    frame_id                                                                      18899
    region_name                                                                  Nevada
    priority                                                                        4.0
    ...
    start_aggressive                                                2001-01-31 00:00:00
    end_aggressive                                                  2001-03-28 00:00:00
    start_median                                                    2001-01-02 00:00:00
    end_median                                                      2001-04-02 00:00:00
    start_conservative                                              2000-11-03 00:00:00
    end_conservative                                                2001-04-17 00:00:00
    blackout_duration_aggressive                                                   56.0
    blackout_duration_median                                                       90.0
    blackout_duration_conservative                                                165.0

    Output format is the following:

        {
        "metadata": {
            "generation_time": "2025-02-13T14:54:35.966369",
            "input_file": "blockout_snow_disp_s1.geojson",
            "output_file": "opera-disp-s1-blackout-dates-2025-02-13.json"
        },
        "blackout_dates": {
            "95": [
            [
                "2016-11-01T00:00:00",
                "2017-05-31T23:59:59"
            ],
            ...
            [
                "2023-11-01T00:00:00",
                "2024-05-31T23:59:59"
            ]
            ],
            "96": [
            [
                "2016-11-01T00:00:00",
                "2017-02-28T23:59:59"
                ...

    Each DISP-S1 frame contains a date range where CSLCs should be excluded
    due to snow cover (or, for the case of Central America, extreme rainy season).

    Parameters
    ----------
    input_file : str
        The path to the input GeoJSON file.
    max_default_duration : float
        The maximum number of days to blackout before switching to the
        more aggressive mode (the mode which keeps more data).
        Default is 180 (6 months).

    Returns
    -------
    dict
        A dictionary containing the blackout dates for each frame_id and metadata.

    Notes
    -----
    This function reads a GeoJSON file, processes it to create blackout dates
    for periods where to_process is 0, saves the result to a JSON file,
    and returns the resulting dictionary.

    """
    generation_time = datetime.now().strftime("%Y-%m-%d")
    output_filename = f"opera-disp-s1-blackout-dates-{generation_time}.json"

    result: dict = {
        "metadata": {
            "generation_time": datetime.now().isoformat(),
            "max_default_duration": max_default_duration,
            "input_file": Path(input_file).name,
            "output_file": output_filename,
        },
    }

    # Load the snow-analysis table
    gdf = gpd.read_parquet(input_file)
    if "start_selected" not in gdf.columns:
        gdf = _select_blackout_dates(gdf, max_default_duration)

    # Span of calendar years to pre-compute
    years = range(2015, 2030)

    blackout_dates: dict[str, list[list[str]]] = {}
    for tup in gdf.itertuples():
        # skip rows missing a valid window
        if pd.isna(tup.start_selected) or pd.isna(tup.end_selected):
            continue

        blackout_dates[str(tup.frame_id)] = _yearly_windows(
            tup.start_selected, tup.end_selected, years
        )

    result["blackout_dates"] = blackout_dates

    with open(output_filename, "w") as f:
        json.dump(result, f, indent=2)

    logger.info(
        "Blackout JSON created: %s (%d frames)", output_filename, len(blackout_dates)
    )
    return result


def _select_blackout_dates(
    gdf: gpd.GeoDataFrame, max_default_duration: float
) -> gpd.GeoDataFrame:
    use_aggressive_mask = gdf.blackout_duration_median > max_default_duration
    gdf["start_selected"] = np.where(
        use_aggressive_mask,
        gdf["start_aggressive"],
        gdf["start_median"],
    )
    gdf["end_selected"] = np.where(
        use_aggressive_mask,
        gdf["end_aggressive"],
        gdf["end_median"],
    )
    gdf["blackout_duration_selected"] = np.where(
        use_aggressive_mask,
        gdf["blackout_duration_aggressive"],
        gdf["blackout_duration_median"],
    )
    gdf["mode_selected"] = np.where(use_aggressive_mask, "aggressive", "median")
    return gdf


def gdf_to_blackout_json(input_file: Path | str) -> dict:
    """Create a JSON of blackout periods for DISP-S1.

    Convert a GeoJSON file with year, month, frame_id, and to_process information
    into a JSON with blackout dates.

    Output format is the following:

        {
        "metadata": {
            "generation_time": "2025-02-13T14:54:35.966369",
            "input_file": "blockout_snow_disp_s1.geojson",
            "output_file": "opera-disp-s1-blackout-dates-2025-02-13.json"
        },
        "blackout_dates": {
            "95": [
            [
                "2016-11-01T00:00:00",
                "2017-05-31T23:59:59"
            ],
            ...
            [
                "2023-11-01T00:00:00",
                "2024-05-31T23:59:59"
            ]
            ],
            "96": [
            [
                "2016-11-01T00:00:00",
                "2017-02-28T23:59:59"
                ...

    Each DISP-S1 frame contains a date range where CSLCs should be excluded
    due to snow cover (or, for the case of Central America, extreme rainy season).

    Parameters
    ----------
    input_file : str
        The path to the input GeoJSON file.

    Returns
    -------
    dict
        A dictionary containing the blackout dates for each frame_id and metadata.

    Notes
    -----
    This function reads a GeoJSON file, processes it to create blackout dates
    for periods where to_process is 0, saves the result to a JSON file,
    and returns the resulting dictionary.

    """
    # Read the GeoJSON file
    gdf = gpd.read_file(input_file)

    blackout_dates: dict[str, list[list[str]]] = {}

    for frame_id, group in gdf.groupby("frame_id"):
        frame_dates = []
        blackout_start = None

        # Sort the group by year and month
        group_sorted = group.sort_values(["year", "month"])

        for _, row in group_sorted.iterrows():
            year = int(row["year"])
            month = int(row["month"])
            to_process = int(row["to_process"])

            current_date = datetime(year, month, 1)

            if to_process == 0:
                if blackout_start is None:
                    blackout_start = current_date
            else:
                if blackout_start is not None:
                    # End the blackout period at the very end of the previous month
                    blackout_end = current_date - timedelta(days=1)
                    blackout_end = blackout_end.replace(hour=23, minute=59, second=59)
                    frame_dates.append(
                        [blackout_start.isoformat(), blackout_end.isoformat()]
                    )
                    blackout_start = None

        # Check if there's an ongoing blackout period at the end of the data
        if blackout_start is not None:
            frame_dates.append(
                [
                    blackout_start.isoformat(),
                    str(group_sorted.iloc[-1]["year"]) + "-12-31",
                ]
            )

        blackout_dates[str(frame_id)] = frame_dates or []

    generation_time = datetime.now().strftime("%Y-%m-%d")
    output_filename = f"opera-disp-s1-blackout-dates-{generation_time}.json"

    result = {
        "metadata": {
            "generation_time": datetime.now().isoformat(),
            "input_file": Path(input_file).name,
            "output_file": output_filename,
        },
        "blackout_dates": blackout_dates,
    }

    # Save the result to a JSON file
    with open(output_filename, "w") as f:
        json.dump(result, f, indent=2)

    logger.info(f"JSON file created: {output_filename}")

    return result


@click.command()
@click.argument("input_file")
@click.option("--max-default-duration", default=240.0)
def create_blackout(
    input_file: Path | str,
    max_default_duration: float = 240.0,
):
    """Create blackout periods JSON for DISP-S1 from `snow_month_filter` outputs."""
    return snow_months_to_blackout_json(
        input_file=input_file,
        max_default_duration=max_default_duration,
    )


if __name__ == "__main__":
    import sys

    try:
        snow_months_to_blackout_json(sys.argv[1])
    except:  # noqa: E722
        gdf_to_blackout_json(sys.argv[1])
