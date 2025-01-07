"""Module to read the OPERA North America shape.

Data comes from:
https://raw.githubusercontent.com/nasa/opera-sds-pcm/refs/heads/develop/geo/opera_NA_expanded.geojson

"""

from pathlib import Path

import geopandas as gpd
from shapely import GeometryType


def get_opera_na_shape() -> GeometryType.MULTIPOLYGON:
    """Read the OPERA North America geometry as a shapely `multipolygon`.

    Data source
    https://raw.githubusercontent.com/nasa/opera-sds-pcm/refs/heads/develop/geo/opera_NA_expanded.geojson

    """
    filename = Path(__file__).parent / "data" / "north_america_opera.geojson.zip"
    na_gpd = gpd.read_file(filename)
    # Combine all geometries in the GeoDataFrame into one MultiPolygon
    return na_gpd.geometry.unary_union
