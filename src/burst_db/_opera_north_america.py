from pathlib import Path

import geopandas as gpd
from shapely import GeometryType


def get_opera_na_shape() -> GeometryType.MULTIPOLYGON:
    filename = Path(__file__).parent / "data" / "north_america_opera.geojson.zip"
    na_gpd = gpd.read_file(filename)
    # Combine all geometries in the GeoDataFrame into one MultiPolygon
    return na_gpd.geometry.unary_union
