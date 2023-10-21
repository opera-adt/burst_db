# Historical burst database creation

Collection of scripts used to generate historical burst catalog.
The three main steps were

1. Use `download_asf_granule_list.py` to query [CMR's STAC catalog](https://cmr.earthdata.nasa.gov/cloudstac/ASF) containing the list of Sentinel-1 granules available on ASF. This creates one list of granules per day.

- Approx runtime: ~5-20 minutes

2. Download the SAFE metadata for each Sentinel-1 url using `download.py`, which relies on [`asfsmd`](https://github.com/avalentino/asfsmd)

- This bulk of the downloads were done on AWS to [enable direct S3 bucket access](https://github.com/scottstanie/asfsmd/tree/s3fs-client), which speeds up the metadata download from ~30 seconds per granule to 1-2 seconds
- Approximately 2 million granules exist on ASF. Runtime depends on ability to scale out API calls.

3. From the Annotation XML files in each granule, run `parse_bursts.py` (which rips out a subset of [s1reader](https://github.com/opera-adt/s1-reader)) to determine which burst IDs are in a granule. Record the burst ID, sensing_time, granule, and approximate WKT geometry

After creating the Burst ID/Frame ID geopackage database and the historical catalog of ~50 million bursts, a merged database was created with `export_denormalized_db.py`.
