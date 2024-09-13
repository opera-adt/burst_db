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

from asfsmd.cli import _get_auth
from eof import download
from tqdm.contrib.concurrent import thread_map

from burst_db.historical_bursts import download_asf_granule_list, parse_bursts
from burst_db.historical_bursts.download_annotations import download_safe_metadata


def _setup_log() -> logging.Logger:
    logger = logging.getLogger("burst_db")
    if logger.hasHandlers():
        logger.handlers.clear()

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
    return logger


logger = _setup_log()


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
    auth = _get_auth()
    logger.info(f"{len(products) = }")
    out_products = [(out_dir / (str(p) + ".SAFE")) for p in products]

    remaining_products = [p for p in out_products if not p.exists()]
    logger.info(f"{len(remaining_products) = }")
    remaining_names = [p.name.replace(".SAFE", "") for p in remaining_products]
    batches = [
        remaining_names[i : i + batch_size]
        for i in range(0, len(remaining_names), batch_size)
    ]

    def _run_download(batch):
        download_safe_metadata(batch, outdir=out_dir, skip_errors=False, auth=auth)

    thread_map(_run_download, batches, max_workers=5)
    return out_products


def get_burst_csvs(
    safe_list: Sequence[str | Path],
    out_dir: Path,
    out_csv: Path,
    orbit_file: Path,
) -> None:
    """Parse the SAFE metadata to get the list of bursts ."""
    logger.info("Finding bursts in SAFE files")
    csv_files = parse_bursts.make_all_safe_metadata(
        safe_list=safe_list,
        out_dir=out_dir,
        orbit_file=orbit_file,
    )
    # Combine all the CSVs into one per date
    parse_bursts._combine_csvs_by_date(csv_files, out_csv, no_clean=False)


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


def main() -> Path:
    """Run the script to generate the list of bursts from a single date."""
    args = _get_cli_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    date_str = args.date.strftime("%Y%m%d")
    logger.info("Getting list of granules for %s from CMR STAC catalog", args.date)
    granule_list_file = get_asf_list(args.date, out_dir)
    granules = granule_list_file.read_text().splitlines()

    out_csv = out_dir / f"{date_str}.csv"
    out_dir_by_date = out_dir / f"{date_str}"
    out_dir_by_date.mkdir(exist_ok=True)

    # Filename for the final
    for satellite in ["A", "B"]:
        if not parse_bursts._is_valid_date(args.date, satellite):
            logger.info("Skipping Sentinel-%s, not valid for %s", satellite, args.date)
            continue

        logger.info("Processing Sentinel-%s", satellite)
        # Create the output directory
        out_dir = args.out_dir

        # Download the SAFE metadata
        annotation_files = get_annotations(granules, out_dir_by_date)

        # Get the orbit for the date
        orbit_file = download.main(
            save_dir=args.orbit_dir,
            date=date_str,
            mission=f"S1{satellite}",
        )[0]

        # Parse the SAFE metadata to get the list of bursts
        cur_burst_csv_file = out_dir_by_date / f"{date_str}_{satellite}.csv"
        get_burst_csvs(
            safe_list=annotation_files,
            out_dir=out_dir_by_date,
            out_csv=cur_burst_csv_file,
            orbit_file=orbit_file,
        )

        # Append this to the total burst CSV output in `out_dir`
        with open(out_csv, "a") as f:
            f.write(cur_burst_csv_file.read_text())

    return out_csv


if __name__ == "__main__":
    main()
