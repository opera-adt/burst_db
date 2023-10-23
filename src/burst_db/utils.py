from __future__ import annotations

import json
import zipfile
from pathlib import Path

from shapely import box


def read_zipped_json(filename: str | Path) -> dict:
    """Read a ".json.zip" file into a dict."""
    with zipfile.ZipFile(filename) as zf:
        bytes = zf.read(str(Path(filename).name).replace(".zip", ""))
        return json.loads(bytes.decode())


def write_zipped_json(json_path: str, dict_out: dict, level: int = 6):
    json_zip_path = str(json_path) + ".zip"
    with zipfile.ZipFile(
        json_zip_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=level
    ) as zf:
        zf.writestr(json_path, json.dumps(dict_out))


def build_wkt_from_bbox(xmin: float, ymin: float, xmax: float, ymax: float) -> str:
    """Convert bounding box coordinates to WKT POLYGON string."""
    return box(xmin, ymin, xmax, ymax).wkt
