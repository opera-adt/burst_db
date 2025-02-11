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


def calculate_reference_dates(
    consistent_json_file: str | None = None,
    desired_month_by_frame: dict[str, int] | None = None,
    interval_years: float = 1.0,
    min_acquisitions_per_batch: int = 15,
) -> dict[str, list[str]]:
    """Generate reference dates for each DISP-S1 Frame.

    If 'desired_month_by_frame' is provided (non-empty), we do a simple
    month-based approach that ignores 'consistent_json_file'. Otherwise,
    we do the original 'once a year from first date' logic reading from
    consistent_json_file.

    Parameters
    ----------
    consistent_json_file : str or None
        Path to the input JSON file with consistent data, or None if we're doing
        month-based references.
    desired_month_by_frame : dict[str, int], optional
        Dictionary of frame IDs -> desired reference month (1..12).
        If given, we do a simple month-based approach.
    interval_years : float, optional
        Approximate interval in years between reference dates (default is 1.0).
    min_acquisitions_per_batch : int, optional
        Minimum number of acquisitions required between reference dates (default is 15).

    Returns
    -------
    dict[str, dict[str, list]]
        Dictionary with frame IDs as keys, containing reference dates
        (and optionally grouped_sensing_times, etc.).

    """
    if desired_month_by_frame:
        # If we have a desired_month_by_frame, do the simple "month-based" approach
        return _generate_month_based_dates(desired_month_by_frame)
    else:
        # Otherwise, do the logic based on the actual data contained in
        # the consistent_json_file
        return _generate_by_consistent(
            consistent_json_file,
            interval_years=interval_years,
            min_acquisitions_per_batch=min_acquisitions_per_batch,
        )


def _generate_month_based_dates(
    desired_month_by_frame: dict[str, int],
    start_year: int = 2016,
    end_year: int = 2030,
) -> dict[str, list[str]]:
    """Generate reference dates on the 1st of the desired month each year.

    Skips additional references for frames in FRAMES_TO_SKIP.
    """
    reference_dates: dict[str, list[str]] = {}

    for frame_id, desired_month in desired_month_by_frame.items():
        all_dates = []
        for year in range(start_year, end_year):
            dt = datetime(year, desired_month, 1)
            dt_str = dt.strftime("%Y-%m-%dT%H:%M:%S")
            all_dates.append(dt_str)

        reference_dates[frame_id] = all_dates

    return reference_dates


def _generate_by_consistent(
    consistent_json_file: str | None = None,
    interval_years: float = 1.0,
    min_acquisitions_per_batch: int = 15,
) -> dict[str, list[str]]:
    if not consistent_json_file:
        raise ValueError(
            "No consistent_json_file provided and no desired_month_by_frame given. "
            "We have no data for references."
        )

    with open(consistent_json_file, "r") as f:
        consistent_data = json.load(f)
        # Check if we have the version with top-level metadata/data
        if "data" in consistent_data:
            consistent_data = consistent_data["data"]

    reference_dates: dict[str, list[str]] = {}
    interval_days = int(interval_years * 365.25)

    for frame_id, frame_data in consistent_data.items():
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
        current_group = []

        for date in sensing_times:
            if not ref_dates:
                ref_dates.append(date)
                current_group = [date]
                continue

            current_group.append(date)

            current_interval = (date - ref_dates[0]).days
            is_interval_passed = current_interval >= len(ref_dates) * interval_days
            is_event_date = date.date() in [ed.date() for ed in frame_event_dates]

            if is_interval_passed or is_event_date:
                if len(current_group) >= min_acquisitions_per_batch:
                    ref_dates.append(date)
                    grouped_sensing_times.append(current_group)
                    current_group = []
                elif is_event_date:
                    # Merge with previous group if not enough acquisitions
                    if grouped_sensing_times:
                        grouped_sensing_times[-1].extend(current_group)
                    ref_dates[-1] = date
                    current_group = []

        # Check for leftover acquisitions at the end
        if current_group:
            grouped_sensing_times.append(current_group)
            if len(current_group) < min_acquisitions_per_batch:
                logger.debug(
                    f"Frame {frame_id} has only {len(current_group)} acquisitions "
                    "in the last batch."
                )

        # Build final structure
        ref_dates_str = [d.strftime("%Y-%m-%dT%H:%M:%S") for d in ref_dates]

        # Account for frames where we skip the reference change
        if str(frame_id) in FRAMES_TO_SKIP:
            # Keep only the first reference date
            ref_dates_str = ref_dates_str[:1]

        reference_dates[frame_id] = ref_dates_str

    return reference_dates


@click.command()
@click.argument("consistent_json_file", required=False, type=click.Path(exists=True))
@click.option(
    "--blackout-file",
    type=click.Path(exists=True),
    default=None,
    help=(
        "Path to blackout JSON file. If given, use month-based references derived"
        " from snow blackouts.",
    ),
)
@click.option(
    "--output",
    help=(
        "Manually name the output JSON file. If not given, defaults to "
        "'opera-disp-s1-reference-dates-{today}-minimal.json'"
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
    consistent_json_file, blackout_file, output, interval, min_acquisitions
):
    r"""Generate a reference dates JSON file for InSAR time series processing.

    \b
    CONSISTENT_JSON_FILE: Path to the input JSON file with consistent data
    (only required if NOT using --blackout-file).
    """
    # Decide if we're doing month-based or the original approach
    desired_month_by_frame = None
    if blackout_file:
        desired_month_by_frame = build_desired_month_map_from_blackout(blackout_file)

    reference_dates = calculate_reference_dates(
        consistent_json_file=consistent_json_file,
        interval_years=interval,
        min_acquisitions_per_batch=min_acquisitions,
        desired_month_by_frame=desired_month_by_frame,
    )

    # Choose output names
    if output is None:
        today_str = datetime.today().strftime("%Y-%m-%d")
        output = f"opera-disp-s1-reference-dates-minimal-{today_str}.json"

    def _dt_to_str(obj):
        return obj.isoformat() if isinstance(obj, datetime) else obj

    total_out = {
        "metadata": {
            "generation_time": datetime.today(),
            "consistent_json_file": consistent_json_file,
            "blackout_file": blackout_file,
            "interval": interval,
            "min_acquisitions": min_acquisitions,
        },
        "data": reference_dates,
    }
    # Write the full output
    with open(output, "w") as f:
        json.dump(total_out, f, indent=2, default=_dt_to_str)
    click.echo(f"Reference dates JSON  created: {output}")


def pick_month_based_on_snow(num_blackouts: int) -> int:
    """Pick a desired month based the number of snow blackout intervals recorded."""
    if num_blackouts == 0:
        return 11  # No snow: November is probably safe
    elif num_blackouts <= 5:
        return 9  # Occasionally sees snow -> pick September
    else:
        return 7  # Lots of snow -> pick July


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
