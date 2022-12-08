#!/usr/bin/env python
'''
Loads ESA burst map in spatialite3; Writeout the attributes for OPERA project to the database

Input: ESA Burst map (for IW in spatialize3 format)
Output: Augmented Burst map (ESA burst map + Burst geogrid information for OPERA)
Optional output: Deployable database (EPSG, xmin, ymin, xmax, ymax only)
'''

import argparse
import os
import shutil
import sqlite3

from osgeo import ogr, osr

import burst_database_core as bd

def get_args():
    '''
    Parse the arguments and return
    '''

    parser = argparse.ArgumentParser(
                description=('Generate burst bounding box DB for OPERA SAS '
                             'from ESA burstmap'),
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('sqlite_path_in',
                        type=str,
                        help='Path to the source ESA burstmap.')

    parser.add_argument('sqlite_path_out',
                        type=str,
                        help='Path to the output sqlite file')

    parser.add_argument('-mxy',
                        type=float,
                        default=[5000, 5000],
                        nargs=2,
                        help='x/y margin of the bounding boxes in meters')

    parser.add_argument('-sxy',
                        type=int,
                        default=[30, 30],
                        nargs=2,
                        help='x/y snapping values in meters for bounding boxes')

    parser.add_argument('-d','--deployable',
                        type=str,
                        help='(Option) Path to the deployable database')

    args = parser.parse_args()

    return args



def main(args):
    '''
    Docstring please
    '''

    if os.path.exists(args.sqlite_path_out):
        print('Output SQLITE file exists. Removing.')
        os.remove(args.sqlite_path_out)

    # Make a copy of the original ESA burst map data; Play with the duplicate.
    shutil.copy(args.sqlite_path_in, args.sqlite_path_out)

    conn = sqlite3.connect(args.sqlite_path_out)
    cur = conn.cursor()

    # Add additional columns to the table
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN burst_id_jpl text')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN EPSG integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmin integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN xmax integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymin integer')
    cur.execute('ALTER TABLE burst_id_map ADD COLUMN ymax integer')

    # Create index onthe table to expedite the search process
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
    print('Determining EPSG for the bursts on IW2')
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
        envelope_geom_tformed = list(geom_burst.GetEnvelope())

        # Apply margin
        envelope_geom_tformed[0] -= args.mxy[0]
        envelope_geom_tformed[1] -= args.mxy[1]
        envelope_geom_tformed[2] += args.mxy[0]
        envelope_geom_tformed[3] += args.mxy[1]

        xmin, ymin, xmax, ymax = bd.snap_extent(tuple(envelope_geom_tformed),
                                                args.sxy[0],
                                                args.sxy[1])

        str_sql = f'UPDATE  burst_id_map SET '\
                  f'xmin={xmin}, xmax={xmax}, ymin={ymin}, ymax={ymax} WHERE '\
                  f'relative_orbit_number={track_in} AND burst_id={burst_id_in_track} AND '\
                  f'subswath_name="{str_subswath}"'

        cur.execute(str_sql)

        # de-reference the objects from feature
        transform = None
        srs_out = None

        # intermediate commit
        if i_row % 1000 == 0:
            conn.commit()

    # Final commit to the database
    conn.commit()
    print('\n')
    del query_result_all, list_all

    # Generate the deployable DB
    if args.deployable:
        if os.path.exists(args.deployable):
            print(f'Deployable DB exists: {args.deployable}. Deleting.')
            os.remove(args.deployable)

        records_burst_data = bd.extract_burst_geogrid_data(args.sqlite_path_out)
        bd.export_to_sqlite(records_burst_data, args.deployable)

    print('\nProcessing completed!')


if __name__=='__main__':
    arg_in = get_args()
    main(arg_in)
