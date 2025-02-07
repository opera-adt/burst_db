import json
import logging
from datetime import datetime

import click

logger = logging.getLogger(__name__)
EVENT_DATES_BY_FRAME = {
    # Ridgecrest: ascending and descending frames
    "16941": ["2019-07-06"],
    "18903": ["2019-07-06"],
}

FRAMES_TO_SKIP = {
    # Hawaii big island:
    "23211",
    "23212",
    "33038",
    "33039",
}


def calculate_reference_dates(
    consistent_json_file: str,
    interval_years: float = 1.0,
    min_acquisitions_per_batch: int = 15,
    desired_month_by_frame: dict[str, int] | None = None,
) -> dict[str, dict[str, list]]:
    """Generate reference dates for each DISP-S1 Frame.

    Parameters
    ----------
    consistent_json_file : str
        Path to the input JSON file with consistent data.
    interval_years : float, optional
        Approximate interval in years between reference dates (default is 1.0).
    min_acquisitions_per_batch : int, optional
        Minimum number of acquisitions required between reference dates (default is 15).
    desired_month_by_frame : dict[str, list[str]], optional
        Dictionary of frame IDs to desired month of reference.


    Returns
    -------
    dict[str, dict[str, list]]
        Dictionary with frame IDs as keys, containing reference dates
        grouped sensing times.

    """
    desired_month_by_frame = desired_month_by_frame or {}
    if desired_month_by_frame:
        return {
            fid: [
                datetime(y, desired_month_by_frame[str(fid)], 1).strftime(
                    "%Y-%m-%dT%H:%M:%S"
                )
                for y in range(2016, 2030)
            ]
            for fid in desired_month_by_frame
        }

    with open(consistent_json_file, "r") as f:
        consistent_data = json.load(f)
        # Check if we have the verison with top level metadata/data
        if "data" in consistent_data:
            consistent_data = consistent_data["data"]

    reference_dates: dict[str, dict[str, list]] = {}
    interval_days = int(interval_years * 365.25)

    for frame_id, frame_data in consistent_data.items():
        desired_month_by_frame.get(str(frame_id), None)
        sensing_times = [
            datetime.strptime(t, "%Y-%m-%dT%H:%M:%S")
            for t in frame_data["sensing_time_list"]
        ]
        frame_event_dates = [
            datetime.strptime(d, "%Y-%m-%d")
            for d in EVENT_DATES_BY_FRAME.get(str(frame_id), [])
        ]

        ref_dates: list[datetime] = []
        grouped_sensing_times = []
        # Tracking how many dates are between the last reference change and current
        # sensing time. gets reset with each reference change
        current_group = []

        for date in sensing_times:
            if not ref_dates:
                ref_dates.append(date)
                current_group = [date]
                continue

            current_group.append(date)

            is_interval_passed = (date - ref_dates[0]).days >= len(
                ref_dates
            ) * interval_days
            is_event_date = date.date() in [ed.date() for ed in frame_event_dates]

            if is_interval_passed or is_event_date:
                if len(current_group) >= min_acquisitions_per_batch:
                    ref_dates.append(date)
                    grouped_sensing_times.append(current_group)
                    current_group = []
                elif is_event_date:
                    # Merge with previous group if not enough acquisitions
                    # TODO: We need to spot check if this leads to a poor reference date
                    # e.g. if there's an "event" in winter, in alaska. We don't want
                    # to set the reference to winter.
                    if grouped_sensing_times:
                        grouped_sensing_times[-1].extend(current_group)
                    ref_dates[-1] = date
                    current_group = []

        # Check for excess at the end
        if current_group:
            grouped_sensing_times.append(current_group)
            if len(current_group) < min_acquisitions_per_batch:
                logger.debug(
                    f"Frame {frame_id} has only {len(current_group)} acquisitions "
                    "in the last batch."
                )

        reference_dates[frame_id] = {
            "reference_dates": [d.strftime("%Y-%m-%dT%H:%M:%S") for d in ref_dates],
            "grouped_sensing_times": [
                [d.strftime("%Y-%m-%dT%H:%M:%S") for d in group]
                for group in grouped_sensing_times
            ],
            "acquisition_counts": [len(group) for group in grouped_sensing_times],
        }

        # Account for frames where we skip the reference change
        if str(frame_id) in FRAMES_TO_SKIP:
            cur_dates = reference_dates[frame_id]["reference_dates"]
            reference_dates[frame_id]["reference_dates"] = cur_dates[:1]

    return reference_dates


@click.command()
@click.argument("consistent_json_file", type=click.Path(exists=True))
@click.option(
    "--output",
    help=(
        "Manually name the output JSON file. If, not given, defaults to"
        " 'opera-disp-s1-reference-dates-{today}.json'"
    ),
)
@click.option(
    "--interval",
    type=float,
    default=1.0,
    help="Approximate nominal interval between reference dates (in years)",
    show_default=True,
)
@click.option(
    "--min-acquisitions",
    type=int,
    default=15,
    help="Minimum number of acquisitions required between reference dates",
    show_default=True,
)
def make_reference_dates(
    consistent_json_file,
    output,
    interval,
    min_acquisitions,
):
    """Generate a reference dates JSON file for InSAR time series processing.

    CONSISTENT_JSON_FILE: Path to the input JSON file with consistent data.
    """
    reference_dates = calculate_reference_dates(
        consistent_json_file, interval, min_acquisitions
    )
    total_out = {
        "metadata": {
            "generation_time": datetime.today(),
            "consistent_json_file": consistent_json_file,
            "interval": interval,
            "min_acquisitions": min_acquisitions,
        },
        "data": reference_dates,
    }
    if output is None:
        today_str = datetime.today().strftime("%Y-%m-%d")
        output_file = f"opera-disp-s1-reference-dates-{today_str}.json"
        minimal_output_file = f"opera-disp-s1-reference-dates-minimal-{today_str}.json"
    else:
        output_file = output
        minimal_output_file = output_file.replace(".json", "-minimal.json")

    def _dt_to_str(obj):
        return obj.isoformat() if isinstance(obj, datetime) else obj

    with open(output_file, "w") as f:
        json.dump(total_out, f, indent=2, default=_dt_to_str)
    click.echo(f"Reference dates JSON with full sensing times created: {output_file}")

    minimal_version = {
        frame_id: data["reference_dates"] for frame_id, data in reference_dates.items()
    }
    minimal_out = total_out.copy()
    minimal_out["data"] = minimal_version
    with open(minimal_output_file, "w") as f:
        json.dump(minimal_out, f, indent=2, default=_dt_to_str)
    click.echo(f"Minimal reference dates JSON file created: {minimal_output_file}")


def pick_month_based_on_snow(num_blackouts: int) -> int:
    """Pick a desired month based the snow blackout intervals recorded."""
    if num_blackouts == 0:
        return 11  # Now snow: November should be safe
    elif num_blackouts <= 5:
        return 9  # Frame occasionally sees snow -> September
    else:
        return 7  # Much snow -> use July


def build_desired_month_map_from_blackout(json_file: str) -> dict[str, int]:
    """Parse the blackouts JSON and decide the best reference month for each frame.

    Parameters
    ----------
    json_file : str
        Path to the JSON file with blackout info.

    Returns
    -------
    dict[str, int]
        Dictionary mapping frame_id (as string) -> desired month (1..12)

    """
    with open(json_file, "r") as f:
        blackout_data = json.load(f)

    # The block of interest is blackout_data["blackout_dates"]
    frames_blackouts = blackout_data["blackout_dates"]

    desired_month_by_frame = {}

    for frame_id_str, intervals in frames_blackouts.items():
        # intervals is a list of [start_time, end_time]
        num_intervals = len(intervals)

        desired_month = pick_month_based_on_snow(num_intervals)
        desired_month_by_frame[frame_id_str] = desired_month

    return desired_month_by_frame
