import json
import logging
from datetime import datetime, timedelta
from pathlib import Path

import geopandas as gpd
import pandas as pd

logger = logging.getLogger(__name__)


def gdf_to_blackout_json(input_file: Path | str) -> dict:
    """Create a JSON of blackout periods for DISP-S1.

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

        # Sort the group by start_date
        group_sorted = group.sort_values("start_date")

        for _, row in group_sorted.iterrows():
            start_date = pd.to_datetime(row["start_date"])
            pd.to_datetime(row["end_date"])
            to_process = int(row["to_process"])

            if to_process == 0:
                if blackout_start is None:
                    blackout_start = start_date
            else:
                if blackout_start is not None:
                    # End the blackout period at the end of the previous day
                    blackout_end = start_date - timedelta(days=1)
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
                    group_sorted.iloc[-1]["end_date"].isoformat(),
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
