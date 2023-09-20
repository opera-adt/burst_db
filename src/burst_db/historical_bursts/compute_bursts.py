#!/usr/bin/env python
from __future__ import annotations

import argparse
import datetime
import os
from dataclasses import dataclass, field
import boto3
import zipfile

# from concurrent.futures import ProcessPoolExecutor
import warnings
from pathlib import Path
from itertools import chain
import shutil

import asf_search as asf
from asfsmd import make_patterns, download_annotations
import backoff
import geopandas as gpd
import pandas as pd
import requests
import s1reader
from joblib import Parallel, delayed
from shapely.ops import unary_union
from shapely.geometry import shape
from tqdm.contrib.concurrent import thread_map
from tqdm.auto import tqdm

import logging
from rich.logging import RichHandler

logger = logging.getLogger("burst_db")
h = RichHandler(rich_tracebacks=True, log_time_format="[%Y-%m-%d %H:%M:%S]")
logger.addHandler(h)
logger.setLevel(logging.INFO)

logging.getLogger("asfsmd").setLevel(logging.WARNING)


def get_burst_rows(
    safe_file: str | Path, out_dir: str | Path, orbit_dir="/home/staniewi/dev/orbits/"
):
    try:
        outfile = (Path(out_dir) / Path(safe_file).stem).with_suffix(".csv")
        if outfile.exists():
            return
        orbit_file = s1reader.get_orbit_file_from_dir(safe_file, orbit_dir=orbit_dir)
        all_rows = []
        for iw in [1, 2, 3]:
            bursts = s1reader.load_bursts(
                safe_file, orbit_file, swath_num=iw, flag_apply_eap=False
            )
            for b in bursts:
                all_rows.append(_to_row(b, safe_file))

        pd.DataFrame(all_rows).to_csv(outfile, header=False, mode="w", index=False)
    except Exception as e:
        print(f"Failure on {safe_file}: {e}")
        outfile = (Path(out_dir) / f"failure_{Path(safe_file).stem}").touch()


def _to_row(burst, safe_file):
    burst_id = burst.burst_id
    dt = burst.sensing_start.isoformat()
    # Get an approximate border, fewer points
    border = unary_union(burst.border).simplify(1 / 3600).wkt
    return burst_id, dt, border, Path(safe_file).stem


def make_all_safe_metadata(
    *,
    out_dir: str | Path,
    dir_with_safes: str | Path | None,
    safe_list: list[str | Path] | None = None,
    orbit_dir="/home/staniewi/dev/orbits/",
    max_jobs=20,
):
    if safe_list is None:
        safe_list = sorted(Path(dir_with_safes).glob("*.SAFE"))
        print(f"Found {len(safe_list)} SAFE dirs in {dir_with_safes}")

    print(f"Writing CSVs to {out_dir}")

    warnings.filterwarnings("ignore", category=UserWarning)  # s1reader is chatty
    Parallel(n_jobs=max_jobs)(
        delayed(get_burst_rows)(f, out_dir=out_dir, orbit_dir=orbit_dir)
        for f in safe_list
    )


def add_to_file(results: asf.ASFProduct, out: Path):
    rdf = gpd.GeoDataFrame.from_features(results.geojson())
    mode = "a" if out.exists() else "w"
    rdf.to_file(out, mode=mode)


@dataclass
class StacSearch:
    start_date: datetime.date = datetime.date(2014, 10, 3)
    end_date: datetime.date = field(default_factory=datetime.date.today)
    max_workers: int = 10
    output_dir: Path = Path(".")

    @backoff.on_exception(backoff.expo, Exception, max_tries=3)
    @staticmethod
    def get_safes_by_date(
        date: datetime.date, missions=["A", "B"], verbose=False
    ) -> list[str]:
        """Get the list of (VV, IW) SAFEs acquired on one date."""
        # sub catalog per date
        date_list_url = "https://cmr.earthdata.nasa.gov/stac/ASF/collections/SENTINEL-1{sat}_SLC.v1/{date_str}"

        safe_names = []
        for sat in missions:
            resp = requests.get(
                date_list_url.format(sat=sat, date_str=date.strftime("%Y/%m/%d"))
            )
            try:
                resp.raise_for_status()
            except requests.HTTPError:
                if verbose:
                    print(f"Failed for {date}. Skipping.")
                continue

            for item in resp.json()["links"]:
                if not item["rel"] == "item":
                    # e.g. 'self', 'root', 'parent'
                    continue
                # Strip the -SLC, which we'll add in later.
                # We want the name to match the normal file names
                s = item["title"].replace("-SLC", "")
                # Check that the beam mode is IW (skip EW, WV, S3)
                if s[4:6] != "IW":
                    continue
                # Only save VV (so SDV or SSV)
                if s[15] != "V":
                    continue
                safe_names.append(s)

        return safe_names

    @backoff.on_exception(backoff.expo, Exception, max_tries=3)
    @staticmethod
    def get_safe_metadata(safe_name: str) -> tuple[str, str]:
        """Get the geojson WKT and concept ID for one SAFE granule."""
        item_url = "https://cmr.earthdata.nasa.gov/stac/ASF/collections/SENTINEL-1{sat}_SLC.v1/items/{safe_name}-SLC"
        # example:
        # https://cmr.earthdata.nasa.gov/stac/ASF/collections/SENTINEL-1A_SLC.v1/items/S1A_IW_SLC__1SDV_20150302T000329_20150302T000356_004845_006086_51B0-SLC
        sat = "A" if safe_name.startswith("S1A") else "B"
        resp = requests.get(item_url.format(safe_name=safe_name, sat=sat))
        resp.raise_for_status()
        js = resp.json()
        # Get the polygon
        polygon = shape(js["geometry"])

        # Get the "Concept" url which has the S3 bucked
        # example: "https://cmr.earthdata.nasa.gov/search/concepts/G1345380785-ASF.json"
        concept_url = [
            link["href"] for link in js["links"] if link["href"].endswith("-ASF.json")
        ][0]
        concept_id = concept_url.split("/")[-1].replace(".json", "")

        return polygon.wkt, concept_id

    def get_all_safe_names(self, overwrite: bool = False, missions=["A", "B"]):
        """Grab the list of SAFE names for each date in parallel."""

        def _save_save_list(date):
            output_name = self.output_dir / f"safes-{date.strftime('%Y-%m-%d')}.txt"
            if output_name.exists() and not overwrite:
                print(f"{output_name} exists, skipping")
            safe_list = StacSearch.get_safes_by_date(date, missions=missions)
            with open(output_name, "a") as f:
                f.write("\n")
                f.write("\n".join(safe_list))
                f.write("\n")

        dates = pd.date_range(self.start_date, self.end_date, freq="1D").to_pydatetime()
        thread_map(_save_save_list, dates, max_workers=self.max_workers)

    def get_all_safe_metadata(self, overwrite: bool = False):
        """After downloading the safe names, query for each SAFE metadata."""

        safe_lists = sorted(self.output_dir.glob("safes-*.txt"))
        print(f"Found {len(safe_lists)} files in {self.output_dir.resolve()}")
        line_counts = {f: _line_count(f) for f in safe_lists}
        print(f"Total number of SAFES: {sum(line_counts.values())}")
        # Divide the outputs into months
        # ...


def _line_count(filename):
    """Fast line count for a text file in python."""

    def _make_gen(reader):
        b = reader(1024 * 1024)
        while b:
            yield b
            b = reader(1024 * 1024)

    with open(filename, "rb") as f:
        f_gen = _make_gen(f.raw.read)
        return sum(buf.count(b"\n") for buf in f_gen)


def download_safe_metadata(
    product_names: list[str],
    pol: str = "vv",
    outdir: str | Path = Path("."),
    skip_if_exists: bool = True,
    use_s3: bool = False,
):
    if use_s3:
        # Set the ASFSMD_CLIENT environment variable to use S3
        os.environ["ASFSMD_CLIENT"] = "s3fs"

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
        logger.error(f"Error downloading data from {remaining_products}", exc_info=True)


@backoff.on_exception(backoff.expo, Exception, max_tries=3)
def _download_safe_metadata(
    product_name: str,
    pol: str = "vv",
    outdir: Path = Path("."),
):
    """Use `asfsmd` to get the SAFE metadata for a product."""

    patterns = make_patterns(pol=pol)
    download_annotations(product_name, patterns=patterns, outdir=outdir)


def zip_and_upload(
    files: list[Path], bucket_name: str, folder_name: str, remove_local: bool = True
):
    """
    Zips a list of files and uploads the resulting zip file to S3.

    Parameters
    ----------
    files : list[Path]
        list of file paths to be zipped.
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

    # Zip the files
    # Create an S3 client
    s3 = boto3.client("s3")
    for file in files:
        zip_file = Path(str(file) + ".zip")
        with zipfile.ZipFile(zip_file, "w") as zipf:
            zipf.write(file, file.name)

        # Upload the zip file to S3
        s3.upload_file(
            Filename=zip_file, Bucket=bucket_name, Key=f"{folder_name}/{zip_file}"
        )

        # Optionally: remove the local zip file after upload
        if remove_local:
            shutil.rmtree(file)
            Path(f"{zip_file}.zip").unlink()


def _get_product_list_cmr(search_dir: str):
    safe_lists = Path(search_dir).glob("safes-*.txt")
    logger.info(f"Found {len(safe_lists)} text files of SAFE products.")
    return list(
        chain.from_iterable(path.read_text().splitlines() for path in safe_lists)
    )


def main() -> None:
    """Download Sentinel-1 metadata from a WKT file."""

    parser = argparse.ArgumentParser(
        description="Download S1 metadata from a list of SAFE granule names.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--out-dir", default="scratch", type=Path, help="Temporary output directory.")
    parser.add_argument("--bucket", help="S3 bucket name to upload results.")
    parser.add_argument("--folder", help="S3 folder name within bucket to upload results.")
    parser.add_argument(
        "--safe-list", help="Text file with list of SAFE granule names.", type=Path
    )
    parser.add_argument(
        "--batch-size",
        default=500,
        type=int,
        help="Number of SAFE granules to download at once before uploading.",
    )

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
    batches = [product_names[i:i + batch_size] for i in range(0, len(product_names), batch_size)]
    logger.info(f"Divided {len(product_names)} products into {len(batches)} batches.")

    # Download the SAFE metadata
    for batch in tqdm(batches):
        download_safe_metadata(batch, outdir=out_dir, use_s3=True)
        cur_files = list(out_dir.glob("*.SAFE"))
        # Upload the files to S3
        if arg_dict["bucket"] and arg_dict["folder"]:
            zip_and_upload(cur_files, arg_dict["bucket"], arg_dict["folder"])



if __name__ == "__main__":
    main()
