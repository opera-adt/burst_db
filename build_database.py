#!/usr/bin/env python
'''
Loads ESA burst map in spatialite3; Writeout the attributes for OPERA project to the database

Input: ESA Burst map (for IW in spatialize3 format)
Output: Augmented Burst map (ESA burst map + Burst geogrid information for OPERA)

TODO: Take another look at the docstring above.
'''
import argparse
import datetime
import json
import os
import shutil
import sqlite3

import numpy as np
from osgeo import ogr, osr

import build_database_core as bd


if __name__=='__main__':
    
    PATH_SRC_ROOT = os.getenv('HOME') + '/Documents/DATA/burst_db'
    PATH_DST_ROOT = os.path.dirname(__file__)

    STR_TIMESTAMP = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    PATH_DATABASE_SRC = f'{PATH_SRC_ROOT}/burst_map_IW_000001_375887.sqlite3'
    PATH_DATABASE_DST = f'{PATH_DST_ROOT}/data/burst_map_IW_000001_375887.OPERA-JPL.'\
                        f'{STR_TIMESTAMP}.sqlite3'

    MARGIN_X = 1000.0
    MARGIN_Y = 1000.0
    SNAP_X = 30.0
    SNAP_Y = 30.0

    # other flag
    export_db_for_opera = True
    # Flag to choose to add bursts' center
    # (approximated from bounding boxes' center) and area of the bounding box
    verbose_database = False

    # Make a copy of the original ESA burst map data; Play with the duplicate.
    shutil.copy(PATH_DATABASE_SRC, PATH_DATABASE_DST)


    conn = sqlite3.connect(PATH_DATABASE_DST)
    cur = conn.cursor()

    # Add additional columns to the table
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN burst_id_jpl text')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN EPSG integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmin integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmax integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymin integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymax integer')

    #cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmin float')
    #cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmax float')
    #cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymin float')
    #cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymax float')

    cur.execute('CREATE INDEX index_burst ON burst_id_map '\
                '(relative_orbit_number, burst_id, subswath_name)')

    # Assign the JPL burst ID on every row
    query_result_all = cur.execute('SELECT * FROM burst_id_map')
    col_id = {}
    for i, column_desc in enumerate(query_result_all.description):
        col_id[column_desc[0]] = i

    list_all = query_result_all.fetchall()
    num_all = len(list_all)

    print('Assigning JPL burst ID to all rows')
    for i_row, row in enumerate(list_all):
        print(f'processing: {i_row + 1:,} / {num_all:,}', end='\r')
        track_in = row[col_id['relative_orbit_number']]
        burst_id_in_track = row[col_id['burst_id']]
        str_subswath = row[col_id['subswath_name']]
        str_burstid_jpl = bd.get_burst_id(track_in, burst_id_in_track, str_subswath)

        str_sql = (f'UPDATE burst_id_map SET burst_id_jpl="{str_burstid_jpl}" WHERE '
                   f'relative_orbit_number={track_in} AND burst_id={burst_id_in_track} AND '
                   f'subswath_name="{str_subswath}"')

        cur.execute(str_sql)

        if i_row % 1000 == 0:
            conn.commit()

    # Final commit to the database
    conn.commit()
    print('\n')
    del query_result_all, list_all


    # Determine EPSG based on the centroid location of the burst in IW2
    print('Determining EPSG for each bursts')
    query_result_iw2 = cur.execute('SELECT * FROM burst_id_map WHERE subswath_name = "IW2"')
    col_id = {}
    for i_col, column_desc in enumerate(query_result_iw2.description):
        col_id[column_desc[0]] = i_col

    list_iw2 = query_result_iw2.fetchall()
    num_iw2 = len(list_iw2)

    for i_iw2, row in enumerate(list_iw2):
        print(f'processing: {i_iw2 + 1:,} / {num_iw2:,}', end=' - ')

        track_in = row[col_id['relative_orbit_number']]
        burst_id_in_track = row[col_id['burst_id']]

        print(f'Track: {track_in}, Burst: {burst_id_in_track}',end='\r')

        geom_burst = ogr.CreateGeometryFromWkb(row[col_id['GEOMETRY']])

        x_centroid, y_centroid = bd.get_centroid_multipolygon(geom_burst)
        epsg_burst = bd.get_point_epsg(y_centroid,
                                       x_centroid)

        str_sql_epsg = (f'UPDATE  burst_id_map SET EPSG={epsg_burst} WHERE '
                        f'relative_orbit_number={track_in} AND burst_id={burst_id_in_track}')

        cur.execute(str_sql_epsg)

        if i_iw2 % 1000 == 0:
            conn.commit()

    conn.commit()
    print('\n Calculating bounding box')
    # Calculate bounding box for every row
    srs_in = osr.SpatialReference()
    srs_in.ImportFromEPSG(4326)
    srs_in.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)


    query_result_all = cur.execute('SELECT * FROM burst_id_map')
    col_id = {}
    for i_col, column_desc in enumerate(query_result_all.description):
        col_id[column_desc[0]] = i_col

    list_all = query_result_all.fetchall()
    num_all = len(list_all)

    for i_row, row in enumerate(list_all):
        print(f'processing: {i_row + 1:,} / {num_all:,}', end='\r')
        track_in = row[col_id['relative_orbit_number']]
        burst_id_in_track = row[col_id['burst_id']]
        str_subswath=row[col_id['subswath_name']]
        epsg_burst=row[col_id['EPSG']]

        geom_burst=ogr.CreateGeometryFromWkb(row[col_id['GEOMETRY']])

        srs_out = osr.SpatialReference()
        srs_out.ImportFromEPSG(epsg_burst)
        transform = osr.CoordinateTransformation(srs_in, srs_out)

        geom_burst.Transform(transform)
        dict_geom_tformed = json.loads(geom_burst.ExportToJson())
        envelope_geom_tformed = geom_burst.GetEnvelope()

        xmin = int(np.round((envelope_geom_tformed[0] - MARGIN_X) / SNAP_X) * SNAP_X)
        xmax = int(np.round((envelope_geom_tformed[1] + MARGIN_X) / SNAP_X) * SNAP_X)
        ymin = int(np.round((envelope_geom_tformed[2] - MARGIN_Y) / SNAP_Y) * SNAP_Y)
        ymax = int(np.round((envelope_geom_tformed[3] + MARGIN_Y) / SNAP_Y) * SNAP_Y)

        str_sql = f'UPDATE  burst_id_map SET '\
                  f'xmin={xmin}, xmax={xmax}, ymin={ymin}, ymax={ymax} WHERE '\
                  f'relative_orbit_number={track_in} AND burst_id={burst_id_in_track} AND '\
                  f'subswath_name="{str_subswath}"'

        cur.execute(str_sql)

        # dereference the objects from feature
        transform = None
        srs_out = None
        nparr_coord_tformed = None

        # intermediate commit
        if i_row % 1000 == 0:
            conn.commit()

    # Final commit to the database
    conn.commit()
    print('\n')
    del query_result_all, list_all

    print('\nProcessing completed!')
