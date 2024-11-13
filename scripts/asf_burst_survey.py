#!/usr/bin/env python
# /// script
# requires-python = ">=3.11"
# dependencies = [
#   "asf_search",
#   "tqdm",
#   "backoff",
# ]
# ///
"""Survey the ASF burst catalog.

Runnable with `pipx`:

pipx run asf_burst_survey.py burst_ids_north_america.csv

To create the list of North america burst ids:

import opera_utils

gdf_bursts = opera_utils.get_burst_id_geojson(True)
burst_ids_unique = (
    gdf_bursts[gdf_bursts.is_north_america == "1"]
    .burst_id_jpl.str.strip("t")
    .str.upper().unique()
)
pd.DataFrame(burst_ids_unique).to_csv("burst_ids_north_america.csv")


"""

from tqdm.contrib.concurrent import thread_map
from itertools import repeat
import argparse
import csv
import backoff
from pathlib import Path
import logging
from typing import List, Optional
from dataclasses import dataclass

import asf_search

# Configure logging
logger = logging.getLogger(__name__)


@dataclass(order=True)
class BurstResult:
    """Container for ASF search results for a burst ID"""

    burst_id: str
    sensing_time: str
    esa_granule: str
    asf_name: str

    def __post_init__(self):
        self.burst_id = self.burst_id.lower()
        if not self.burst_id.startswith("t"):
            self.burst_id = "t" + self.burst_id


@backoff.on_exception(backoff.expo, Exception, max_tries=5)
def search_burst(burst_id: str, max_results: int = 2000) -> List[BurstResult]:
    """
    Search ASF for a given burst ID.

    Parameters
    ----------
    burst_id : str
        The burst ID to search for
    max_results : int, optional
        Maximum number of results to return

    Returns
    -------
    List[BurstResult]
        List of search results
    """
    results = asf_search.search(
        dataset="SLC-BURST", fullBurstID=burst_id, maxResults=max_results
    )

    return [
        BurstResult(
            burst_id=burst_id,
            sensing_time=result.properties["burst"]["azimuthTime"],
            esa_granule=result.umm["InputGranules"][0],
            asf_name=result.properties["fileID"],
        )
        for result in results
    ]


def save_results(results: List[BurstResult], output_path: Path):
    """Save results to CSV file"""
    if not results:
        return

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(
            f, fieldnames=["burst_id", "sensing_time", "esa_granule", "asf_name"]
        )
        writer.writeheader()
        for result in sorted(results):
            writer.writerow(vars(result))


def process_burst(
    burst_id: str, output_dir: Path, max_results: int
) -> Optional[List[BurstResult]]:
    """Process a single burst ID with progress tracking"""
    try:
        results = search_burst(burst_id, max_results)
        if results:
            save_results(results, output_dir / f"{burst_id.replace('/', '_')}.csv")
        return results
    except Exception as e:
        logger.error(f"Error processing burst {burst_id}: {e}")
        return None


def main():
    parser = argparse.ArgumentParser(
        description="Search ASF for burst IDs and save results",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input_file", type=Path, help="File containing burst IDs (one per line)"
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("burst_results"),
        help="Directory to store output CSV files",
    )

    parser.add_argument(
        "--max-results", type=int, default=2000, help="Maximum results per burst ID"
    )

    parser.add_argument(
        "--workers", type=int, default=5, help="Number of concurrent workers"
    )

    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Enable verbose logging"
    )

    args = parser.parse_args()

    # Configure logging
    logger.setLevel(logging.INFO)
    logger.addHandler(logging.StreamHandler())

    # Create output directory
    args.output_dir.mkdir(exist_ok=True, parents=True)

    # Read burst IDs
    with open(args.input_file) as f:
        burst_ids = [
            line.strip().upper().replace("T", "") for line in f if line.strip()
        ]

    logger.info(f"Processing {len(burst_ids)} burst IDs with {args.workers} workers")

    results = thread_map(
        process_burst,
        burst_ids,
        repeat(args.output_dir),
        repeat(args.max_results),
        max_workers=args.workers,
    )

    logger.info(
        f"Processing complete. Found {len(results)} results across {len(burst_ids)} burst IDs"
    )


if __name__ == "__main__":
    main()
