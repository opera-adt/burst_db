import json
import logging
from datetime import datetime
from pathlib import Path

import click
import opera_utils.download
from opera_utils import group_by_date

from burst_db.utils import batched

logger = logging.getLogger("burst_db")


OUTPUT_TYPES = ["granule", "https", "s3"]


def get_urls_for_frame(
    frame_id: str | int,
    json_file: str,
    output_type: str = OUTPUT_TYPES[0],
    start_date: datetime | None = None,
    end_date: datetime | None = None,
):
    """Retrieve URLs for a DISP-S1 Frame ID using the provided JSON file."""
    if output_type not in OUTPUT_TYPES:
        raise ValueError(
            f"Unknown output_type: {output_type}, must be in {OUTPUT_TYPES}"
        )

    data: dict[str, dict[str, list[str]]]
    with open(json_file, "r") as f:
        loaded = json.load(f)
        # Account for the "unnested" version, or nested withing data/metadata keys
        data = loaded.get("data", loaded)

    if str(frame_id) not in data:
        raise ValueError(f"Frame {frame_id} not found in the JSON file.")

    frame_data = data[str(frame_id)]
    burst_ids = frame_data.get("burst_id_list", [])
    sensing_times = frame_data.get("sensing_time_list", [])
    # Convert sensing times to dates and filter based on start_date and end_date
    sensing_date_set = {
        datetime.fromisoformat(t).date()
        for t in sensing_times
        if (start_date is None or datetime.fromisoformat(t).date() >= start_date.date())
        and (end_date is None or datetime.fromisoformat(t).date() <= end_date.date())
    }
    expected_total_files = len(burst_ids) * len(sensing_date_set)

    logger.info(f"Searching {json_file} for bursts for Frame {frame_id}.")
    if not sensing_times:
        raise ValueError(f"No sensing times found for Frame {frame_id}.")
    else:
        logger.info(f"Found {len(sensing_times)} sensing times for Frame {frame_id}.")

    if not burst_ids:
        raise ValueError(f"No burst IDs found for Frame {frame_id}.")
    else:
        logger.info(f"Found {len(burst_ids)} burst IDs for Frame {frame_id}:")
        logger.info(",".join(burst_ids))

    logger.info("Performing ASF search for the specified burst IDs.")
    try:
        search_results = opera_utils.download.search_cslcs(
            burst_ids=burst_ids, start=start_date, end=end_date
        )
    except Exception as e:
        raise ValueError("ASF search failed") from e

    if isinstance(search_results, tuple):
        search_results, _ = search_results  # Ignore missing data options

    logger.info(f"Total search results before filtering: {len(search_results)}")

    # Filter results by sensing time
    logger.info("Filtering results by sensing time.")
    filtered_results = []
    for result in search_results:
        burst_id = result.properties.get("operaBurstID", "").lower()
        start_time_str = result.properties.get("startTime", "")
        try:
            start_time = datetime.fromisoformat(start_time_str)
            date_str = start_time.date()
        except ValueError:
            logger.warning(
                f"Invalid startTime format for burst ID {burst_id}: {start_time_str}"
            )
            continue

        if date_str in sensing_date_set:
            filtered_results.append(result)
        else:
            logger.debug(f"Skipping {burst_id}: {date_str}")

    logger.info(
        f"Search results after filtering by sensing time: {len(filtered_results)}"
    )

    if not filtered_results:
        logger.warning("No search results match the provided sensing times.")
        return
    # Apply deduplication
    deduped_results = opera_utils.download.filter_results_by_date_and_version(
        filtered_results
    )
    # logger.info(f"Results after deduplication: {len(deduped_results)}")
    if len(deduped_results) != expected_total_files:
        result_str = "\n".join(
            [r.properties.get("fileName", "") for r in deduped_results]
        )
        logger.error(f"Unexpected results: {result_str}")
        raise ValueError(
            f"Expected to find {expected_total_files}, found {len(deduped_results)}"
        )

    # Depending on output_type, extract the desired information
    if output_type.lower() == "granule":
        return [r.properties.get("fileName", "") for r in deduped_results]

    return opera_utils.download.get_urls(deduped_results, type_=output_type.lower())


@click.command()
@click.argument("frame_id", type=str)
@click.option(
    "--json-file",
    type=click.Path(exists=True, dir_okay=False, readable=True),
    default="opera-disp-s1-consistent-burst-ids-2016-07-01_to_2024-09-04.json",
    show_default=True,
    help="Path to the JSON file containing frame ID mappings.",
)
@click.option(
    "--output-type",
    type=click.Choice(["https", "granule", "s3"], case_sensitive=False),
    default="https",
    show_default=True,
    help="Type of output to display.",
)
@click.option(
    "--start-date",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    help="Start date for filtering results (YYYY-MM-DD).",
)
@click.option(
    "--end-date",
    type=click.DateTime(formats=["%Y-%m-%d"]),
    help="End date for filtering results (YYYY-MM-DD).",
)
@click.option(
    "--ministack-size",
    type=int,
    default=None,
    show_default=True,
    help="Number of URLs per ministack.",
)
@click.option(
    "--output-dir",
    type=click.Path(file_okay=False, writable=True),
    default=".",
    show_default=True,
    help="Directory to save ministack files.",
)
def urls_for_frame(
    frame_id: str,
    json_file: str,
    output_type: str,
    start_date: datetime | None,
    end_date: datetime | None,
    ministack_size: int | None,
    output_dir: str,
):
    """Retrieve URLs for a specific FRAME_ID using the provided JSON file.

    FRAME_ID: The frame ID to retrieve URLs for (e.g., "831").
    """
    urls = sorted(get_urls_for_frame(
        frame_id=frame_id,
        json_file=json_file,
        output_type=output_type,
        start_date=start_date if start_date else None,
        end_date=end_date if end_date else None,
    ))
    if not ministack_size:
        ministack_size = len(urls)

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    date_to_urls = group_by_date(urls, date_idx=0)
    for i, batch in enumerate(batched(date_to_urls.items(), ministack_size)):
        output_file = output_path / f"ministack-{i:02d}.txt"
        num_urls = 0
        with output_file.open("w") as f:
            for _date_tup, cur_urls in batch:
                num_urls += len(cur_urls)
                f.write("\n".join(cur_urls) + "\n")
        click.echo(f"Written {num_urls} URLs to {output_file}")
