#!/usr/bin/env python
"""Functions for downloading the list of Sentinel-1 granules available from ASF.

Uses the STAC catalog run by CMR at https://cmr.earthdata.nasa.gov/cloudstac/ASF/
"""
from __future__ import annotations

import argparse
import datetime
import logging
from dataclasses import dataclass, field
from pathlib import Path

import backoff
import pandas as pd
import requests
from shapely.geometry import shape
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
        date_list_url = "https://cmr.earthdata.nasa.gov/cloudstac/ASF/collections/SENTINEL-1{sat}_SLC.v1/{date_str}"  # noqa

        safe_names = []
        for sat in missions:
            resp = requests.get(
                date_list_url.format(sat=sat, date_str=date.strftime("%Y/%m/%d"))
            )
            try:
                resp.raise_for_status()
            except requests.HTTPError:
                if verbose:
                    logger.warning(f"Failed for {date}. Skipping.")
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
        item_url = "https://cmr.earthdata.nasa.gov/stac/ASF/collections/SENTINEL-1{sat}_SLC.v1/items/{safe_name}-SLC"  # noqa
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

    def get_all_safe_names(
        self, overwrite: bool = False, missions=["A", "B"]
    ) -> list[Path]:
        """Grab the list of SAFE names for each date in parallel."""

        def _save_save_list(date) -> Path:
            output_name = self.output_dir / f"safes-{date.strftime('%Y-%m-%d')}.txt"
            if output_name.exists() and not overwrite:
                logger.info(f"{output_name} exists, skipping")
                return output_name
            safe_list = StacSearch.get_safes_by_date(date, missions=missions)
            with open(output_name, "a") as f:
                f.write("\n")
                f.write("\n".join(safe_list))
                f.write("\n")
            return output_name

        dates = pd.date_range(self.start_date, self.end_date, freq="1D").to_pydatetime()
        if len(dates) > 1:
            return thread_map(_save_save_list, dates, max_workers=self.max_workers)
        else:
            return [_save_save_list(dates[0])]


def main() -> None:
    """Download Sentinel-1 metadata from a WKT file."""

    parser = argparse.ArgumentParser(
        description="Download the available Sentinel-1 SAFE granules from ASF.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--out-dir", default="scratch", type=Path, help="Temporary output directory."
    )
    parser.add_argument(
        "--start-date",
        default=datetime.date(2014, 10, 3),
        type=datetime.date.fromisoformat,
        help="Start date for search.",
    )
    parser.add_argument(
        "--end-date",
        default=datetime.date.today(),
        type=datetime.date.fromisoformat,
        help="End date for search.",
    )
    parser.add_argument(
        "--max-workers",
        default=10,
        type=int,
        help="Number of workers to use for downloading SAFE metadata.",
    )

    args = parser.parse_args()

    # Create the output directory
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Set up the StacSearch object
    stac_search = StacSearch(
        start_date=args.start_date,
        end_date=args.end_date,
        max_workers=args.max_workers,
        output_dir=out_dir,
    )

    # Get the list of SAFE granule names
    stac_search.get_all_safe_names()


if __name__ == "__main__":
    main()
