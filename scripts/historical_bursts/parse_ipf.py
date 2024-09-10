#!/usr/bin/env python
import argparse
import zipfile
from glob import glob
from itertools import islice  # , batched
from pathlib import Path
from typing import Sequence

from tqdm.contrib.concurrent import thread_map


def get_ipf_from_zip(filename: str | Path):
    with zipfile.ZipFile(filename) as zf:
        # zobj = zf.open(f'{p.name.replace(".zip", "")}/manifest.safe')
        # Find manifest.safe in namelist
        for name in zf.namelist():
            if "manifest.safe" in name:
                break
        else:
            raise ValueError(f"Could not find manifest.safe in {filename}")
        zobj = zf.open(name)
        for line_bytes in zobj:
            line = line_bytes.decode().strip()
            if line.startswith('<safe:software name="Sentinel-1 IPF"'):
                return line.split("version=")[1].strip('"/>')


def get_ipf_from_zip_list(zip_list: Sequence[str | Path], max_workers: int = 4):
    file_gen = (z for z in zip_list if Path(z).is_file())
    num_files = len(zip_list)
    return thread_map(
        get_ipf_from_zip,
        file_gen,
        max_workers=max_workers,
        desc="Getting IPF versions",
        total=num_files,
    )


def batched(iterable, n):
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def save_granule_ipfs(
    zip_list: Sequence[str | Path],
    output_file: Path,
    max_workers: int = 4,
    batch_size: int = 1000,
):
    # Write the header
    output_file.write_text("granule,ipf\n")
    # Iterate over the batches
    for sub_list in batched(zip_list, batch_size):
        # Get the IPF versions
        ipf_list = get_ipf_from_zip_list(sub_list, max_workers=max_workers)
        # Write the results
        for filename, ipf in zip(sub_list, ipf_list):
            granule = Path(filename).name.replace(".zip", "").replace(".SAFE", "")
            with output_file.open("a") as f:
                f.write(f"{granule},{ipf}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("zip_list", nargs="+", help="List of zip files to parse")
    parser.add_argument(
        "-o",
        "--output_file",
        type=Path,
        default="granule_ipf_versions.csv",
        help="Output file",
    )
    parser.add_argument(
        "--max_workers", type=int, default=4, help="Number of workers to use"
    )
    args = parser.parse_args()

    if len(args.zip_list) == 1 and "*" in args.zip_list[0]:
        zip_list = glob(args.zip_list[0])
    else:
        zip_list = args.zip_list

    save_granule_ipfs(
        zip_list, output_file=args.output_file, max_workers=args.max_workers
    )


if __name__ == "__main__":
    main()
