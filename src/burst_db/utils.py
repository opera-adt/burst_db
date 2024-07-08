from __future__ import annotations

import json
import logging
import time
import zipfile
from functools import wraps
from pathlib import Path
from typing import Callable, ParamSpec, TypeVar

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


# Used for callable types
T = TypeVar("T")
P = ParamSpec("P")


def log_runtime(f: Callable[P, T]) -> Callable[P, T]:
    """Decorate a function to time how long it takes to run.

    Usage
    -----
    @log_runtime
    def test_func():
        return 2 + 4
    """
    logger = logging.getLogger(__name__)

    @wraps(f)
    def wrapper(*args: P.args, **kwargs: P.kwargs):
        t1 = time.time()

        result = f(*args, **kwargs)

        t2 = time.time()
        elapsed_seconds = t2 - t1
        elapsed_minutes = elapsed_seconds / 60.0
        time_string = (
            f"Total elapsed time for {f.__module__}.{f.__name__} : "
            f"{elapsed_minutes:.2f} minutes ({elapsed_seconds:.2f} seconds)"
        )

        logger.info(time_string)

        return result

    return wrapper
