"""Module to read the OPERA North America shape.

Data comes from:
https://github.com/OPERA-Cal-Val/DSWx-Validation-Experiments/blob/387acb91d7cce3ec62fac47394c219c74b9dee82/marshak/north_america_boundaries/north_america_opera_with_sub_zone.geojson

"""

from pathlib import Path

import geopandas as gpd
from shapely import GeometryType


def get_opera_na_shape() -> GeometryType.MULTIPOLYGON:
    """Read the OPERA North America geometry as a shapely `multipolygon`.

    Data source
    https://github.com/OPERA-Cal-Val/DSWx-Validation-Experiments/blob/387acb91d7cce3ec62fac47394c219c74b9dee82/marshak/north_america_boundaries/north_america_opera_with_sub_zone.geojson
    """
    filename = Path(__file__).parent / "data" / "north_america_opera.geojson.zip"
    na_gpd = gpd.read_file(filename)
    # Combine all geometries in the GeoDataFrame into one MultiPolygon
    return na_gpd.geometry.unary_union
