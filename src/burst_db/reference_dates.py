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
    consistent_json_file: str,
    interval_years: float = 1.0,
    min_acquisitions_per_batch: int = 15,
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
    event_dates_by_frame : dict[str, list[str]], optional
        Dictionary of frame IDs to lists of event dates to include as reference dates.

    Returns
    -------
    dict[str, dict[str, list]]
        Dictionary with frame IDs as keys, containing reference dates
        grouped sensing times.

    """
    with open(consistent_json_file, "r") as f:
        consistent_data = json.load(f)

    reference_dates: dict[str, dict[str, list]] = {}
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

    return reference_dates


@click.command()
@click.argument("consistent_json_file", type=click.Path(exists=True))
@click.option(
    "--output",
    help=(
        "Manually name the output JSON file. If, not given, defaults to"
        " 'opera-disp-s1-consistent-burst-ids-{today}.json'"
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
    OUTPUT_FILE: Path to the output JSON file.
    """
    reference_dates = calculate_reference_dates(
        consistent_json_file, interval, min_acquisitions
    )
    if output is None:
        today_str = datetime.today().strftime("%Y-%m-%d")
        output_file = f"opera-disp-s1-reference-dates-{today_str}.json"
        minimal_output_file = f"opera-disp-s1-reference-dates-minimal-{today_str}.json"
    else:
        output_file = output
        minimal_output_file = output_file.replace(".json", "-minimal.json")

    with open(output_file, "w") as f:
        json.dump(reference_dates, f, indent=2)
    click.echo(f"Reference dates JSON with full sensing times created: {output_file}")

    minimal_version = {
        frame_id: data["reference_dates"] for frame_id, data in reference_dates.items()
    }
    with open(minimal_output_file, "w") as f:
        json.dump(minimal_version, f, indent=2)
    click.echo(f"Minimal reference dates JSON file created: {minimal_output_file}")
