#!/usr/bin/env python
"""Script to generate the list of bursts from a single date.

Required packages are listed in `conda-env.yml` in this directory.
"""
from __future__ import annotations

import argparse
import datetime
import logging
from pathlib import Path
from typing import Sequence

from eof import download
from tqdm.contrib.concurrent import thread_map

from burst_db.historical_bursts import download_asf_granule_list, parse_bursts
from burst_db.historical_bursts.download_annotations import download_safe_metadata

logger = logging.getLogger("burst_db")
# Make a logger good for AWS logs
h = logging.StreamHandler()
h.setFormatter(
    logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
    )
)
logger.addHandler(h)
logger.setLevel(logging.INFO)


def get_asf_list(
    query_date: datetime.date,
    out_dir: Path,
) -> Path:
    """Download the list of available Sentinel-1 granules from ASF.

    Parameters
    ----------
    start_date : datetime.date
        Start date for search.
    end_date : datetime.date
        End date for search.
    out_dir : Path
        Output directory for the list of granules.

    Returns
    -------
    Path
        Path to the list of granules for the given date.
    """
    # Create the output directory
    out_dir.mkdir(parents=True, exist_ok=True)
    # Set up the StacSearch object
    stac_search = download_asf_granule_list.StacSearch(
        start_date=query_date,
        end_date=query_date,
        max_workers=1,
        output_dir=out_dir,
    )

    # Get the list of SAFE granule names
    return stac_search.get_all_safe_names()[0]


def get_annotations(
    products: Sequence[str | Path], out_dir: Path, batch_size: int = 50
):
    """Download the XML annotation files for a list of SAFE products."""
    print(f"{len(products) = }")
    out_products = [Path(out_dir) / p for p in products]

    remaining_products = [
        p for p in products if not (out_dir / (str(p) + ".SAFE")).exists()
    ]
    print(f"{len(remaining_products) = }")
    batches = [
        remaining_products[i : i + batch_size]
        for i in range(0, len(remaining_products), batch_size)
    ]

    def _run_download(batch):
        download_safe_metadata(batch, outdir=out_dir, skip_errors=False)

    thread_map(_run_download, batches, max_workers=5)
    return out_products


def get_burst_csvs(
    date: datetime.date,
    safe_list: Sequence[str | Path],
    out_dir: Path,
    orbit_file: Path,
) -> Path:
    """Parse the SAFE metadata to get the list of bursts ."""
    logger.info("Finding bursts in SAFE files")
    csv_files = parse_bursts.make_all_safe_metadata(
        safe_list=safe_list,
        out_dir=out_dir,
        orbit_file=orbit_file,
    )

    # Combine all the CSVs into one per date
    return parse_bursts._combine_csvs_by_date(csv_files, date, out_dir, no_clean=False)


def _get_cli_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download the available Sentinel-1 SAFE granules from ASF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "-d",
        "--date",
        type=datetime.date.fromisoformat,
        help="Start date for search.",
        required=True,
    )
    parser.add_argument(
        "--out-dir", default="scratch", type=Path, help="Temporary output directory."
    )
    parser.add_argument(
        "--orbit-dir",
        default="scratch",
        type=Path,
        help=(
            "Directory containing existing orbit file. If not provided, will download"
            " precise orbits for `date`."
        ),
    )
    parser.add_argument(
        "--max-workers",
        default=10,
        type=int,
        help="Number of workers to use for downloading SAFE metadata.",
    )

    return parser.parse_args()


def main() -> list[Path]:
    """Run the script to generate the list of bursts from a single date."""

    args = _get_cli_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    granule_list_file = get_asf_list(args.date, out_dir)
    granules = granule_list_file.read_text().splitlines()

    out_paths = []
    for satellite in ["A", "B"]:
        logger.info("Processing Sentinel-%s", satellite)
        # Create the output directory
        out_dir = args.out_dir
        cur_dir = out_dir / f"S1{satellite}"
        cur_dir.mkdir(exist_ok=True)

        # Download the SAFE metadata
        annotation_dir = cur_dir / "annotations"
        annotation_dir.mkdir(exist_ok=True)
        annotation_files = get_annotations(granules, annotation_dir)

        # Get the orbit for the date
        orbit_file = download.main(
            save_dir=args.orbit_dir,
            date=args.date,
            mission=f"S1{satellite}",
        )[0]

        # Parse the SAFE metadata to get the list of bursts
        burst_csv_file = get_burst_csvs(
            date=args.date.date(),
            safe_list=annotation_files,
            out_dir=cur_dir,
            orbit_file=orbit_file,
        )
        out_paths.append(burst_csv_file)
    return out_paths


if __name__ == "__main__":
    main()
