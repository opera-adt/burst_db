from pathlib import Path

import geopandas as gpd
import opera_utils
import pandas as pd

if __name__ == "__main__":
    _here = Path(__file__).parent
    df_volcanoes = pd.read_csv(_here / "USGS_Volcano_Table_2__subset_1-57_.csv")
    gpd.GeoDataFrame(
        df_volcanoes,
        geometry=gpd.points_from_xy(df_volcanoes.longitude, df_volcanoes.latitude),
    )
    gdf_frames = opera_utils.get_frame_geojson(as_geodataframe=True)
    gdf_volcanoes = gpd.GeoDataFrame(
        df_volcanoes,
        geometry=gpd.points_from_xy(df_volcanoes.longitude, df_volcanoes.latitude),
        crs="epsg:4326",
    )
    gdf_priority = gpd.read_file(
        _here / "NApriorityrollout_framebased_v8_13Mar2025.geojson"
    )
    gdf = pd.merge(
        gdf_frames,
        gdf_priority[["frame_id", "priority"]],
        left_on="frame_id",
        right_on="frame_id",
    ).set_crs("epsg:4326")
    volcano_frames = (
        gdf[gdf.priority == 4]
        .sjoin(gdf_volcanoes, how="inner", predicate="intersects")
        .frame_id.unique()
    )
