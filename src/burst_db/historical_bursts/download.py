#!/usr/bin/env python
from __future__ import annotations

import argparse
import logging
import os
import shutil
import time
from itertools import chain
from pathlib import Path

import boto3
from asfsmd import download_annotations, make_patterns

# from eof import products
from tqdm.contrib.concurrent import thread_map

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

logging.getLogger("asfsmd").setLevel(logging.WARNING)


def download_safe_metadata(
    product_names: list[str],
    pol: str = "vv",
    outdir: str | Path = Path("."),
    skip_if_exists: bool = True,
):
    # out_product = (Path(outdir) / product_name).with_suffix(".SAFE")
    out_products = [Path(outdir) / p for p in product_names]
    # keep only the ones that don't exist
    if skip_if_exists:
        remaining_products = [
            p for p, out in zip(product_names, out_products) if not out.exists()
        ]
    else:
        remaining_products = product_names

    logger.info(f"Downloading {len(remaining_products)} products")
    try:
        _download_safe_metadata(remaining_products, pol=pol, outdir=Path(outdir))
    except Exception:
        logger.error(f"Error downloading data from {remaining_products}")
        # Try once more
        time.sleep(20)
        try:
            _download_safe_metadata(remaining_products, pol=pol, outdir=Path(outdir))
        except Exception:
            logger.error(
                f"Error downloading data from {remaining_products}", exc_info=True
            )


# @backoff.on_exception(backoff.expo, Exception, max_tries=2)
def _download_safe_metadata(
    product_names: list[str],
    pol: str = "vv",
    outdir: Path = Path("."),
):
    """Use `asfsmd` to get the SAFE metadata for a product."""

    patterns = make_patterns(pol=pol)
    download_annotations(product_names, patterns=patterns, outdir=outdir)


def zip_and_upload(
    safe_dirs: list[Path], bucket_name: str, folder_name: str, remove_local: bool = True
):
    """
    Zips a list of safe_dirs and uploads the resulting zip file to S3.

    Parameters
    ----------
    safe_dirs : list[Path]
        list of SAFE directory paths to be zipped.
    bucket_name : str
        Name of the S3 bucket to upload to.
    folder_name : str
        S3 folder to upload the zip file to.
    remove_local : bool
        Delete the local zip file and safe folder after upload.

    Returns
    -------
    None
    """

    # Zip the safe_dirs
    # Create an S3 client
    s3 = boto3.client("s3")
    logger.info(
        f"Uploading {len(safe_dirs)} safe_dirs to s3://{bucket_name}/{folder_name}"
    )
    for safe_dir in safe_dirs:
        zip_file = Path(
            shutil.make_archive(
                str(safe_dir),
                "zip",
                root_dir=str(safe_dir.parent),
                base_dir=safe_dir.name,
            )
        )

        # Get the time to use for the S3 key
        key = f"{folder_name}/{zip_file.name}"
        # # Use the year/month like YYYY/MM for a folder structure
        # acq_time = products.Sentinel(str(zip_file)).start_time
        # date_str = acq_time.strftime("%Y/%m")
        # key = f"{folder_name}/{date_str}/{zip_file.name}"

        # Upload the zip safe_dir to S3
        s3.upload_file(Filename=zip_file, Bucket=bucket_name, Key=key)

        # Optionally: remove the local zip safe_dir after upload
        if remove_local:
            logger.debug(f"Removing {safe_dir}")
            shutil.rmtree(safe_dir)
            zip_file.unlink()


def _get_product_list_cmr(search_dir: str):
    safe_lists = list(Path(search_dir).glob("safes-*.txt"))
    logger.info(f"Found {len(safe_lists)} text files of SAFE products.")
    return list(
        chain.from_iterable(path.read_text().splitlines() for path in safe_lists)
    )


def _get_parser():
    parser = argparse.ArgumentParser(
        description="Download S1 metadata from a list of SAFE granule names.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--out-dir", default="scratch", type=Path, help="Temporary output directory."
    )
    parser.add_argument("--bucket", help="S3 bucket name to upload results.")
    parser.add_argument(
        "--folder", help="S3 folder name within bucket to upload results."
    )
    parser.add_argument(
        "--safe-list", help="Text file with list of SAFE granule names.", type=Path
    )
    parser.add_argument(
        "--batch-size",
        default=100,
        type=int,
        help="Number of SAFE granules to download at once before uploading.",
    )
    parser.add_argument(
        "--start-idx",
        default=0,
        type=int,
        help="Index of the first batch to download.",
    )
    parser.add_argument(
        "--end-idx",
        type=int,
        help="Index of the last batch to download.",
    )

    parser.add_argument(
        "--max-workers",
        default=3,
        type=int,
        help=(
            "Number of workers to use for downloading SAFE metadata. "
            "Each worker will get on of the batches of SAFE granules."
        ),
    )
    return parser


def main() -> None:
    """Download Sentinel-1 metadata directly from S3."""
    parser = _get_parser()
    args = parser.parse_args()
    arg_dict = vars(args)

    # Create the output directory
    out_dir = arg_dict["out_dir"]
    out_dir.mkdir(parents=True, exist_ok=True)

    # Get the list of SAFE granule names
    if arg_dict["safe_list"]:
        product_names = arg_dict["safe_list"].read_text().splitlines()
    else:
        product_names = _get_product_list_cmr(out_dir)

    # Divide the list into batches
    batch_size = arg_dict["batch_size"]
    batches = [
        product_names[i : i + batch_size]
        for i in range(0, len(product_names), batch_size)
    ]
    logger.info(f"Divided {len(product_names)} products into {len(batches)} batches.")

    # if arg_dict["use_s3"]:
    # Set the ASFSMD_CLIENT environment variable to use S3
    os.environ["ASFSMD_CLIENT"] = "s3fs"

    def _download_and_save_s3(batch_idx: int):
        logger.info(f"Downloading batch {batch_idx} of {len(batches)}")
        batch = batches[batch_idx]
        cur_out_dir = out_dir / f"batch_{batch_idx}"
        cur_out_dir.mkdir(exist_ok=True)
        download_safe_metadata(batch, outdir=cur_out_dir)
        cur_files = list(cur_out_dir.glob("*.SAFE"))
        # Upload the files to S3
        zip_and_upload(cur_files, arg_dict["bucket"], arg_dict["folder"])

    start_idx = args.start_idx
    end_idx = len(batches) if args.end_idx is None else args.end_idx
    end_idx = min(end_idx, len(batches))

    # Download the SAFE metadata
    thread_map(
        _download_and_save_s3,
        range(start_idx, end_idx),
        desc="Downloading SAFE metadata",
        total=end_idx - start_idx,
        max_workers=arg_dict["max_workers"],
    )
    # )
    # for batch_idx in range(len(batches)):
    #     _download_and_save_s3(batch_idx)


if __name__ == "__main__":
    main()

# screen
# micromamba activate -n base
# total_lines=$(wc -l < remaining_safes.txt)
# lines_per_file=$((total_lines / 5))
# split -l $lines_per_file remaining_safes.txt -d remaining_safes_part_
# ASFSMD_CLIENT=s3fs python download.py --batch 100 --max-work 20 \
#       --bucket burst-database-scott --folder safe_folders_zipped2 \
#       --safe-list ./remaining_safes_part_04
