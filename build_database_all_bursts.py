"""
Loads ESA burst map in spatialite3; Writeout the attributes for OPERA project to the database

Input: ESA Burst map (for IW in spatialize3 format)
Output: Augmented Burst map (ESA burst map + Burst geogrid information for OPERA)
"""

import datetime
import json
import os
import shutil
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import geopandas as gpd

from osgeo import ogr, osr
from shapely import wkb

import build_database as bd




if __name__ == "__main__":
    PATH_SRC_ROOT = sys.argv[1]
    PATH_DST_ROOT = os.path.dirname(__file__)

    STR_TIMESTAMP = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")

    PATH_DATABASE_SRC = f"{PATH_SRC_ROOT}/sqlite/burst_map_IW_000001_375887.sqlite3"
    PATH_DATABASE_DST = (
        f"{PATH_DST_ROOT}/data/burst_map_IW_000001_375887.OPERA-JPL."
        f"{STR_TIMESTAMP}.sqlite3"
    )
    # conn = sqlite3.connect(PATH_DATABASE_DST)

    MARGIN_X = 1000.0
    MARGIN_Y = 1000.0
    SNAP_X = 50.0
    SNAP_Y = 50.0

    TESTING = False

    conn = sqlite3.connect(PATH_DATABASE_SRC)
    print("Loading all rows")
    query = "SELECT * FROM burst_id_map"
    if TESTING:
        query += " LIMIT 1000"
    df = pd.read_sql_query(query, conn)

    geo_series = gpd.GeoSeries.from_wkb(df.GEOMETRY)
    df = gpd.GeoDataFrame(df, geometry=geo_series, crs=4326)

    print("Forming JPL ids")
    burst_id_jpl = (
        "t"
        + df["relative_orbit_number"].astype(str).str.zfill(3)
        + "_"
        + df["burst_id"].astype(str).str.zfill(6)
        + "_"
        + df["subswath_name"]
    )
    df["burst_id_jpl"] = burst_id_jpl

    # Now iterate over the bursts and get the geogrid information
    print("Finding EPSGs and bboxes")
    epsg_bbox = df[['geometry']].apply(bd.get_epsg_bbox, axis=1, result_type="expand") 
    df[['EPSG', 'xmin', 'ymin', 'xmax', 'ymax']] = epsg_bbox

#     print("Loading geometries with OGR")
#     multipolys = df["GEOMETRY"].apply(ogr.CreateGeometryFromWkb)

#     print("Finding centroids")
#     centroids = multipolys.apply(bd.get_centroid_multipolygon)
#     print("Finding EPSG codes from centroids")
#     df["EPSG"] = centroids.apply(
#         lambda lonlat: bd.get_point_epsg(lonlat[1], lonlat[0])
#     )
# # df[['xmin', 'ymin', 'xmax', 'ymax']] = df[['multipolys']].apply(lambda x: x.multipolys.bounds, axis=1, result_type="expand")
#     # list(transform(transform(geom_idl, 32660).centroid, 4326, 32660).coords)
#     # Out[246]: [(179.75286533347008, 52.14289367674507)]

#     # Get the bounding boxes in UTM
#     print("Loading geometries with shapely")
#     df["geom_shapely"] = df["GEOMETRY"].apply(wkb.loads)
#     print("Getting UTM bounding boxes")
#     bboxes = df[["geom_shapely", "EPSG"]].apply(
#         lambda row: bd.get_utm_bbox(
#             row.geom_shapely, row.EPSG, snap=SNAP_X, margin=MARGIN_X
#         ),
#         axis=1,
#     )
#     bboxes = np.array(bboxes.to_list())
#     df["xmin"] = bboxes[:, 0]
#     df["ymin"] = bboxes[:, 1]
#     df["xmax"] = bboxes[:, 2]
#     df["ymax"] = bboxes[:, 3]
    
    df = df.rename({'GEOMETRY': 'geometry_wkb'}, axis=1)
    df = df.drop(['geometry', 'OGC_FID', 'burst_id', 'subswath_name'], axis=1)
    print("Writing to database")
    bd.make_spatialite_db(df, table_name="burst_id_map", outname=PATH_DATABASE_DST)

    # # Make a copy of the original ESA burst map data; Play with the duplicate.
    # Path(f"{PATH_DST_ROOT}/data/").mkdir(exist_ok=True)
    # # shutil.copy(PATH_DATABASE_SRC, PATH_DATABASE_DST)
    # cur = conn.cursor()

    # # # Add additional columns to the table
    # # cur.execute('ALTER TABLE burst_id_map ADD COLUMN burst_id_jpl text')
    # # cur.execute('ALTER TABLE burst_id_map ADD COLUMN EPSG integer')
    # # cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmin float')
    # # cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmax float')
    # # cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymin float')
    # # cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymax float')

    # # cur.execute('CREATE INDEX index_burst ON burst_id_map '\
    # #             '(relative_orbit_number, burst_id, subswath_name)')

    # # Assign the JPL burst ID on every row
    # query_result_all = cur.execute("SELECT * FROM burst_id_map")
    # col_id = {}
    # for i, column_desc in enumerate(query_result_all.description):
    #     col_id[column_desc[0]] = i

    # list_all = query_result_all.fetchall()
    # num_all = len(list_all)

    # for i_row, row in enumerate(list_all):
    #     print(f"processing: {i_row + 1:,} / {num_all:,}", end="\r")
    #     track_in = row[col_id["relative_orbit_number"]]
    #     burst_id_in_track = row[col_id["burst_id"]]
    #     str_subswath = row[col_id["subswath_name"]]
    #     str_burstid_jpl = bd.get_burst_id(track_in, burst_id_in_track, str_subswath)

    #     str_sql = (
    #         f'UPDATE burst_id_map SET burst_id_jpl="{str_burstid_jpl}" WHERE '
    #         f"relative_orbit_number={track_in} AND burst_id={burst_id_in_track} AND "
    #         f'subswath_name="{str_subswath}"'
    #     )

    #     cur.execute(str_sql)

    #     if i_row % 1000 == 0:
    #         conn.commit()

    # # Final commit to the database
    # conn.commit()
    # print("\n")
    # del query_result_all, list_all

    # # Determine EPSG based on the centroid location of the burst in IW2
    # print("Determining EPSG for each bursts")
    # query_result_iw2 = cur.execute(
    #     'SELECT * FROM burst_id_map WHERE subswath_name = "IW2"'
    # )
    # col_id = {}
    # for i_col, column_desc in enumerate(query_result_iw2.description):
    #     col_id[column_desc[0]] = i_col

    # list_iw2 = query_result_iw2.fetchall()
    # num_iw2 = len(list_iw2)

    # for i_iw2, row in enumerate(list_iw2):
    #     print(f"processing: {i_iw2 + 1:,} / {num_iw2:,}", end=" - ")

    #     track_in = row[col_id["relative_orbit_number"]]
    #     burst_id_in_track = row[col_id["burst_id"]]

    #     print(f"Track: {track_in}, Burst: {burst_id_in_track}", end="\r")

    #     geom_burst = ogr.CreateGeometryFromWkb(row[col_id["GEOMETRY"]])

    #     x_centroid, y_centroid = bd.get_centroid_multipolygon(geom_burst)
    #     epsg_burst = bd.get_point_epsg(y_centroid, x_centroid)

    #     str_sql_epsg = (
    #         f"UPDATE  burst_id_map SET EPSG={epsg_burst} WHERE "
    #         f"relative_orbit_number={track_in} AND burst_id={burst_id_in_track}"
    #     )

    #     cur.execute(str_sql_epsg)

    #     if i_iw2 % 1000 == 0:
    #         conn.commit()

    # conn.commit()

    # print("\n Calculating bounding box")
    # # Calculate bounding box for every row
    # srs_in = osr.SpatialReference()
    # srs_in.ImportFromEPSG(4326)
    # srs_in.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)

    # query_result_all = cur.execute("SELECT * FROM burst_id_map")
    # col_id = {}
    # for i_col, column_desc in enumerate(query_result_all.description):
    #     col_id[column_desc[0]] = i_col

    # list_all = query_result_all.fetchall()
    # num_all = len(list_all)

    # for i_row, row in enumerate(list_all):
    #     print(f"processing: {i_row + 1:,} / {num_all:,}", end="\r")
    #     track_in = row[col_id["relative_orbit_number"]]
    #     burst_id_in_track = row[col_id["burst_id"]]
    #     str_subswath = row[col_id["subswath_name"]]
    #     epsg_burst = row[col_id["EPSG"]]

    #     geom_burst = ogr.CreateGeometryFromWkb(row[col_id["GEOMETRY"]])

    #     srs_out = osr.SpatialReference()
    #     srs_out.ImportFromEPSG(epsg_burst)
    #     transform = osr.CoordinateTransformation(srs_in, srs_out)

    #     geom_burst.Transform(transform)
    #     dict_geom_tformed = json.loads(geom_burst.ExportToJson())
    #     # TODO: this is wrong, because it's grabbing the first polygon
    #     # the bounding box will be the bounds of both of them
    #     nparr_coord_tformed = np.array(dict_geom_tformed["coordinates"][0][0])

    #     xmin = np.round((nparr_coord_tformed[:, 0].min() - MARGIN_X) / SNAP_X) * SNAP_X
    #     xmax = np.round((nparr_coord_tformed[:, 0].max() + MARGIN_X) / SNAP_X) * SNAP_X
    #     ymin = np.round((nparr_coord_tformed[:, 1].min() - MARGIN_Y) / SNAP_Y) * SNAP_Y
    #     ymax = np.round((nparr_coord_tformed[:, 1].max() + MARGIN_Y) / SNAP_Y) * SNAP_Y

    #     str_sql = (
    #         f"UPDATE  burst_id_map SET "
    #         f"xmin={xmin}, xmax={xmax}, ymin={ymin}, ymax={ymax} WHERE "
    #         f"relative_orbit_number={track_in} AND burst_id={burst_id_in_track} AND "
    #         f'subswath_name="{str_subswath}"'
    #     )

    #     cur.execute(str_sql)

    #     # dereference the objects from feature
    #     transform = None
    #     srs_out = None
    #     nparr_coord_tformed = None

    #     # intermediate commit
    #     if i_row % 1000 == 0:
    #         conn.commit()

    # # Final commit to the database
    # conn.commit()
    # # Create an index on the string id column
    # cur.execute("CREATE INDEX idx_burst_id_jpl on burst_id_map (burst_id_jpl);")
    # conn.commit()
    # print("\n")
    # del query_result_all, list_all

    # print("\nProcessing completed!")
