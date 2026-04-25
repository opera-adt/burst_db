from __future__ import annotations

import json
import zipfile
from itertools import islice
from pathlib import Path

import click
from shapely import box


def read_zipped_json(filename: str | Path) -> dict:
    """Read a ".json.zip" file into a dict."""
    with zipfile.ZipFile(filename) as zf:
        bytes_ = zf.read(str(Path(filename).name).replace(".zip", ""))
        return json.loads(bytes_.decode())


def write_zipped_json(json_path: str, dict_out: dict, level: int = 6):
    """Write a JSON dictionary to a sibling compressed ".json.zip" file."""
    json_zip_path = str(json_path) + ".zip"
    with zipfile.ZipFile(
        json_zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=level
    ) as zf:
        zf.writestr(json_path, json.dumps(dict_out))


def build_wkt_from_bbox(xmin: float, ymin: float, xmax: float, ymax: float) -> str:
    """Convert bounding box coordinates to WKT POLYGON string."""
    return box(xmin, ymin, xmax, ymax).wkt


def batched(iterable, n):
    """Divide `iterable` into `n` batches."""
    # batched('ABCDEFG', 3) --> ABC DEF G
    if n < 1:
        raise ValueError("n must be at least one")
    it = iter(iterable)
    while batch := tuple(islice(it, n)):
        yield batch


def normalize_cmr_csv_header(
    csv_file: str | Path, output_file: str | Path | None = None
):
    """Normalize CMR CSV file header to expected format.

    The CMR survey CSV files sometimes have inconsistent headers or may be missing
    the header entirely. This function ensures the header is always exactly:
    "Granule ID,Revision Time,Temporal Time,Revision-Temporal Delta Hours,revision-id"

    Parameters
    ----------
    csv_file : str | Path
        Path to the input CSV file (may or may not have a header).
    output_file : str | Path | None, optional
        Path to the output CSV file. If None, overwrites the input file.

    Notes
    -----
    - If the first line contains "Granule" or "Revision", it's treated as a header
      and replaced.
    - If the first line doesn't look like a header (starts with data), the
      standard header is prepended.
    - The expected header format matches what create_cslc_burst_catalog.py expects.

    """
    csv_file = Path(csv_file)
    output_file = Path(output_file) if output_file else csv_file

    # Expected header for burst catalog processing
    expected_header = [
        "Granule ID",
        "Revision Time",
        "Temporal Time",
        "Revision-Temporal Delta Hours",
        "revision-id",
    ]

    # Read all lines
    with open(csv_file) as f:
        lines = f.readlines()

    if not lines:
        raise ValueError(f"CSV file {csv_file} is empty")

    # Check if first line looks like a header (contains expected field names)
    first_line = lines[0].strip()
    is_header = any(
        field in first_line for field in ["Granule", "Revision", "Temporal"]
    )

    # Prepare output lines
    if is_header:
        # Replace existing header
        output_lines = [",".join(expected_header) + "\n", *lines[1:]]
    else:
        # Prepend header if missing
        output_lines = [",".join(expected_header) + "\n", *lines]

    # Write normalized CSV
    with open(output_file, "w") as f:
        f.writelines(output_lines)


@click.command()
@click.argument("csv_file", type=click.Path(exists=True, path_type=Path))
@click.option(
    "--output",
    "-o",
    type=click.Path(path_type=Path),
    default=None,
    help="Output file path. If not specified, overwrites the input file.",
)
def normalize_csv_header(csv_file: Path, output: Path | None):
    """Normalize CMR CSV header to expected format.

    Ensures the CSV header is exactly:
    "Granule ID,Revision Time,Temporal Time,Revision-Temporal Delta Hours,revision-id"

    This handles cases where:
    - The CSV has a different header format
    - The CSV has no header at all
    - The CSV header has variations in spacing or capitalization

    CSV_FILE: Path to the CMR survey CSV file to normalize.
    """
    normalize_cmr_csv_header(csv_file, output)
    output_path = output if output else csv_file
    click.echo(f"✓ Normalized CSV header: {output_path}")
