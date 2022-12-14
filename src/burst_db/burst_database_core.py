'''
A code to build the database for S1 burst coverage
The codes in the `__main__` namespace will be deprecated.
'''

import datetime
import json
import os
from collections import OrderedDict

import sqlite3
import numpy as np
from osgeo import ogr, osr

from burst_db import __version__ as burst_database_version

def get_point_epsg(lat, lon):
    '''
    Get EPSG code based on latitude and longitude
    coordinates of a point
    Copied from geogrid.py in OPERA RTC

    Parameters
    ----------
    lat: float
        Latitude coordinate of the point
    lon: float
        Longitude coordinate of the point

    Returns
    -------
    epsg: int
        UTM zone
    '''

    # "wrap" the longitude value into range [-180.0, 180.0]
    if (lon >= 180.0) or (lon <= -180.0):
        lon = (lon + 180.0) % 360.0 - 180.0

    if lat >= 75.0:
        return 3413
    elif lat <= -60.0:
        return 3031
    elif lat > 0:
        return 32601 + int(np.round((lon + 177) / 6.0))
    elif lat < 0:
        return 32701 + int(np.round((lon + 177) / 6.0))
    else:
        err_str = "'Could not determine EPSG for {0}, {1}'.format(lon, lat))"
        raise ValueError(err_str)


def get_list_polygon_wkt(path_shp: str, epsg_out: str) -> dict:
    '''Take the path to the .shp file.
       Returns the list of the polygons in the input as WKT string
    '''

    drv_in = ogr.GetDriverByName("ESRI Shapefile")
    datasrc_in = drv_in.Open(path_shp, 0)
    lyr_in = datasrc_in.GetLayer()
    srs_in = lyr_in.GetSpatialRef()

    srs_out = osr.SpatialReference()
    srs_out.ImportFromEPSG(epsg_out)
    transform = osr.CoordinateTransformation(srs_in, srs_out)

    num_feat = lyr_in.GetFeatureCount()
    num_field = lyr_in.GetLayerDefn().GetFieldCount()

    dict_out = {}
    dict_out['wkt'] = [None] * num_feat
    dict_out['field'] = [None] * num_field

    #Load the field info and set up the list for the attribute
    for i_field in range(num_field):
        str_field = lyr_in.GetLayerDefn().GetFieldDefn(i_field).GetName()
        dict_out['field'][i_field] = str_field
        dict_out[str_field] = [None]*num_feat

    for i_feat, feat in enumerate(lyr_in):
        print(f'{i_feat + 1} / {num_feat}', end='\r')
        geom = feat.GetGeometryRef()
        geom.Transform(transform)

        dict_out['wkt'][i_feat] = geom.ExportToWkt()
        for str_field in dict_out['field']:
            dict_out[str_field][i_feat] = feat.GetField(str_field)

    print('\n')

    return dict_out


def snap_extent(mimnax_xy: tuple, snap_x: int, snap_y: int):
    '''
    Snap the coordinates in the tuple

    Parameters:
    -----------
    minmax_xy: tuple
        The corner coordinates of the bounding box
        (xmin, ymin, xmax, ymax)
    snap_x: int
        Snap value in x coordinates
    snap_y: int
        Snap value in y coordinates

    Return:
    -------
    _: tuple
        Snapped value of `minmax_xy`

    '''
    # inside tuple: (xmin, ymin, xmax, ymax)
    # i.e. Same as -te option in gdalwarp

    if snap_x > 0:
        xmin_snap = int(np.floor(mimnax_xy[0] / snap_x) * snap_x)
        xmax_snap = int(np.ceil(mimnax_xy[2] / snap_x) * snap_x)
    else:
        xmin_snap = mimnax_xy[0]
        xmax_snap = mimnax_xy[2]

    if snap_y > 0:
        ymin_snap = int(np.floor(mimnax_xy[1] / snap_y) * snap_y)
        ymax_snap = int(np.ceil(mimnax_xy[3] / snap_y) * snap_y)
    else:
        ymin_snap = mimnax_xy[1]
        ymax_snap = mimnax_xy[3]

    return (xmin_snap, ymin_snap, xmax_snap, ymax_snap)


def wkt2extent(str_wkt: str,
               margin_x: float, margin_y: float,
               snap_x: int=0, snap_y: int=0):
    '''
    Calculate the extend of the input polygon, with margin applied.
    Performs snapping when snap>0

    Parameters:
    -----------
    str_wkt: str
        Polygon as WKT string

    margin_x, margin_y: float
        Margins in x / y coordinates to ba added to the polygon's extent [m]

    snap: int
        Snap interval. The x/y coordinates will be
        rounded to the multiples to this value.


    Return:
    -------
    extent: list
        list of float numbers for bounding box
        [xmin, ymin, xmax, ymax]

    '''

    #print(str_wkt)
    token_coord = str_wkt.split('((')[-1]\
                         .split('))')[0]\
                         .split(',')

    num_node = len(token_coord)

    arr_x = np.zeros(num_node)
    arr_y = np.zeros(num_node)
    arr_z = np.zeros(num_node)

    #Test the first coord string to deterimine its dimension
    dimension_coord = len(token_coord[0].split(' '))

    for i_node, str_coord in enumerate(token_coord):
        coord_float = [float(coord) for coord in str_coord.split(' ')]
        arr_x[i_node] = coord_float[0]
        arr_y[i_node] = coord_float[1]

        if dimension_coord > 2:
            arr_z[i_node] = coord_float[2]

    xmin_coord = arr_x.min() - margin_x
    xmax_coord = arr_x.max() + margin_x
    ymin_coord = arr_y.min() - margin_y
    ymax_coord = arr_y.max() + margin_y

    extent = snap_extent((xmin_coord, ymin_coord,
                          xmax_coord, ymax_coord),
                          snap_x, snap_y)

    return extent


def get_burst_id(track:int, burst:int, swath:str) -> str:
    '''Get the string of burst ID.
    '''
    form_burst_id = 't{TRACK:03d}_{BURST:06d}_{SWATH}'
    str_burst_id = form_burst_id.format(TRACK=track,
                                        BURST=burst,
                                        SWATH=swath.lower())

    return str_burst_id


def generate_shp_out(path_shp_in: str, path_shp_out: str,
                     margin_x:float=0.0, margin_y:float=0.0,
                     snap_x:float=5.0, snap_y:float=5.0):
    '''
    Generate a Shapefile whose feature has burst polygon as geometry,
    bounding box information as well as other relevant attributes

    Parameters:
    -----------
    path_shp_in: str
        Subset of burst map in .shp file
    path_shp_out: str
        Output .shp file
    margin_x: float
        Margin to be added in x axis to determine the bounding box [m]
    margin_y: float
        Margin to be added in y axis to determine the bounding box [m]
    snap_x: float
        Snap value to which the x coordinates will be rounded [m]
    snap_y: float
        Snap value to which the y coordinates will be rounded [m]

    '''
    drv_in = ogr.GetDriverByName("ESRI Shapefile")
    datasrc_in = drv_in.Open(path_shp_in, 0)
    lyr_in = datasrc_in.GetLayer()
    srs_in = lyr_in.GetSpatialRef()
    num_feat = lyr_in.GetFeatureCount()

    #set up the output .shp file
    drv_out = ogr.GetDriverByName("ESRI Shapefile")
    datasrc_out = drv_out.CreateDataSource(path_shp_out)
    lyr_out = datasrc_out.CreateLayer('Bursts', srs_in, ogr.wkbPolygon)

    field_burst_id = ogr.FieldDefn("burst_id", ogr.OFTString)
    field_burst_id.SetWidth(22)
    lyr_out.CreateField(field_burst_id)
    lyr_out.CreateField(ogr.FieldDefn('Track', ogr.OFTInteger))
    lyr_out.CreateField(ogr.FieldDefn('EPSG', ogr.OFTInteger))
    lyr_out.CreateField(ogr.FieldDefn('xmin', ogr.OFTReal))
    lyr_out.CreateField(ogr.FieldDefn('xmax', ogr.OFTReal))
    lyr_out.CreateField(ogr.FieldDefn('ymin', ogr.OFTReal))
    lyr_out.CreateField(ogr.FieldDefn('ymax', ogr.OFTReal))

    for i_feat, feat_in in enumerate(lyr_in):
        print(f'Processing: {i_feat + 1} / {num_feat}', end='\r')
        geom = feat_in.GetGeometryRef()
        str_dict_centroid = geom.Centroid().ExportToJson()
        dict_centroid = json.loads(str_dict_centroid)
        epsg_burst = get_point_epsg(dict_centroid['coordinates'][1],
                                    dict_centroid['coordinates'][0])

        track = feat_in.GetField('relative_o')
        butst_id_in_track = feat_in.GetField('burst_id')
        swath = feat_in.GetField('subswath_n')

        str_burst_id = get_burst_id(track, butst_id_in_track, swath)
        srs_out = osr.SpatialReference()
        srs_out.ImportFromEPSG(epsg_burst)
        transform = osr.CoordinateTransformation(srs_in, srs_out)

        wkt_polygon_before_transform = geom.ExportToWkt()

        # Extract the coordinates after the transformation
        geom.Transform(transform)
        dict_geom_tformed = json.loads(geom.ExportToJson())
        nparr_coord_tformed = np.array(dict_geom_tformed['coordinates'][0])

        tuple_extent = ((nparr_coord_tformed[:,0].min() - margin_x),
                        (nparr_coord_tformed[:,0].max() + margin_x),
                        (nparr_coord_tformed[:,1].min() - margin_y),
                        (nparr_coord_tformed[:,1].max() + margin_y))
        xmin, xmax, ymin, ymax = snap_extent(tuple_extent, snap_x, snap_y)

        feat_out = ogr.Feature(lyr_out.GetLayerDefn())
        feat_out.SetField('burst_id', str_burst_id)
        feat_out.SetField('Track', track)
        feat_out.SetField('EPSG', epsg_burst)
        feat_out.SetField('xmin', xmin)
        feat_out.SetField('xmax', xmax)
        feat_out.SetField('ymin', ymin)
        feat_out.SetField('ymax', ymax)
        feat_out.SetGeometry(
            ogr.CreateGeometryFromWkt(wkt_polygon_before_transform))

        lyr_out.CreateFeature(feat_out)

        feat_out = None
        transform = None

    datasrc_out = None


def get_centroid_multipolygon(geometry_in):
    '''
    Calculate the centroid of multipolygon.
    Takes care of geometries separated on +/- 180 degree longitude line

    Parameter:
    ----------
    geometry_in: osgeo.ogr.Geometry
        Burst multipolygon as OSGEO geometry object

    Return:
    -------
    x_centroid, y_centroid: float
        x / y coordinates of the centroid

    '''

    # Detect the number of polygons in the geometry
    dict_geometry = json.loads(geometry_in.ExportToJson())
    num_polygon = len(dict_geometry['coordinates'])

    if num_polygon == 0:
        raise ValueError('Cannot find polygons in the input geometry.')

    elif num_polygon == 1:
        geometry_centroid = geometry_in.Centroid()
        dict_centroid = json.loads(geometry_centroid.ExportToJson())
        x_centroid = dict_centroid['coordinates'][0]
        y_centroid = dict_centroid['coordinates'][1]

    elif num_polygon > 1:
        offset_circular = 360.0
        xy_weight_centroid = np.zeros((num_polygon, 3))

        for id_polygon, nodes_polygon in enumerate(dict_geometry['coordinates']):
            dict_sub_polygon={
                "type": "MultiPolygon",
                "coordinates":[nodes_polygon]}

            json_sub_polygon = json.dumps(dict_sub_polygon)
            geom_sub_polygon = ogr.CreateGeometryFromJson(json_sub_polygon)

            centroid_sub_polygon = geom_sub_polygon.Centroid()
            dict_centroid_sub_polygon = json.loads(
                                            centroid_sub_polygon.ExportToJson())
            x_centroid = dict_centroid_sub_polygon['coordinates'][0]
            y_centroid = dict_centroid_sub_polygon['coordinates'][1]
            area_sub_polygon = geom_sub_polygon.Area()

            xy_weight_centroid[id_polygon,:] = [
                (x_centroid + offset_circular) % offset_circular,
                y_centroid,
                area_sub_polygon]

        # Weighted sum
        x_centroid_weighted_raw = \
            np.sum(xy_weight_centroid[:,0] * xy_weight_centroid[:,2])\
                   / np.sum(xy_weight_centroid[:,2])
        y_centroid_weighted_raw = \
            np.sum(xy_weight_centroid[:,1] * xy_weight_centroid[:,2])\
                   / np.sum(xy_weight_centroid[:,2])

        # Refine the raw result of the weighted sum coordinates
        if x_centroid_weighted_raw > 180.0:
            x_centroid = x_centroid_weighted_raw - offset_circular
        elif x_centroid_weighted_raw < -180.0:
            x_centroid = x_centroid_weighted_raw + offset_circular
        else:
            x_centroid = x_centroid_weighted_raw

        if y_centroid_weighted_raw > 180.0:
            y_centroid = y_centroid_weighted_raw - offset_circular
        elif y_centroid_weighted_raw < -180.0:
            y_centroid = y_centroid_weighted_raw + offset_circular
        else:
            y_centroid = y_centroid_weighted_raw

    return x_centroid, y_centroid


def extract_burst_geogrid_data(path_augmented_burst_map: str):
    '''
    Extract the burst geogrid data from augmented burst map

    Parameters:
    -----------
    path_augmented_burst_map: str
        Path to the augmented burt map

    Return:
    records_out: list
        Burst geogrid information

    '''

    conn = sqlite3.connect(path_augmented_burst_map)
    curs = conn.cursor()
    curs.execute('SELECT burst_id_jpl, EPSG, xmin, ymin, xmax, ymax '
                 'FROM burst_id_map')
    records_out = curs.fetchall()

    curs.close()
    conn.close()

    return records_out


def export_to_csv(records_out: list, path_csv_out: str):
    '''
    Write out the burst geogrid information as .csv file

    Parameters:
        records_out: list
            List of burst geogrid information retrieved from
            `extract_burst_geogrid_data()`
        path_csv_out: str
            path to the .csv file to write out

    '''
    num_record = len(records_out)
    with open(path_csv_out,'w+',encoding='utf-8') as fout:
        # write CSV header
        fout.write('burst_id, EPSG, xmim, ymin, xmax, ymax\n')

        for i_record, record in enumerate(records_out):
            print(f'Writing: {i_record+1} / {num_record}', end='\r')
            fout.write(f'{record[0]}, {record[1]}, {record[2]}, ',
                    f'{record[3]}, {record[4]}, {record[5]}\n' )
        print('\n')


def export_to_json(records_out: list, path_json_out: str):
    '''
    Write out the burst geogrid information as .json file

    Parameters:
        records_out: list
            List of burst geogrid information retrieved from
            `extract_burst_geogrid_data()`
        path_json_out: str
            path to the .json file to write out

    '''

    num_record = len(records_out)

    # Convert the data into dict
    dict_export={}
    for i_record in range(num_record):
        dict_export[records_out[i_record][0]] = {
            'EPSG':records_out[i_record][1],
            'extent':[records_out[i_record][2], records_out[i_record][3],
                    records_out[i_record][4], records_out[i_record][5]]
        }

    with open(path_json_out, 'w+', encoding='utf8') as fout:
        json.dump(dict_export, fout, indent=2)


def export_to_sqlite(records_out: list, path_sqlite_out: str, create_index: bool=False):
    '''
    Write out the burst geogrid information as sqlite database file

    Parameters:
        records_out: list
            List of burst geogrid information retrieved from
            `extract_burst_geogrid_data()`
        path_sqlite_out: str
            path to the .sqlite file to write out

    '''

    num_record = len(records_out)

    # Export to SQLITE
    print('Exporting to SQLITE')
    with sqlite3.connect(path_sqlite_out) as conn_out:
        curs_out = conn_out.cursor()
        curs_out.execute('CREATE TABLE IF NOT EXISTS burst_id_map ('
                         'burst_id_jpl text PRIMARY KEY, EPSG integer, '
                         'xmin integer, ymin integer, xmax integer, ymax integer);')

        for i_record, record in enumerate(records_out):
            print(f' Processing: {i_record:,} / {num_record:,}      ', end='\r')
            str_sql_command = ( 'INSERT INTO burst_id_map '
                                '(burst_id_jpl ,EPSG, xmin, ymin, xmax, ymax) '
                               f'VALUES("{record[0]}", {record[1]}, '
                               f'{record[2]}, {record[3]}, '
                               f'{record[4]}, {record[5]})')

            curs_out.execute(str_sql_command)
        print('\n')
        if create_index:
            curs_out.execute('CREATE INDEX index_burst ON burst_id_map (burst_id_jpl)')

        conn_out.commit()


def writeout_metadata(sql_path, args):
    '''
    Write the metadata (setting, version, creation date, etc.) into db file

    Parameters:
    -----------
    sql_path: str
        Path to the SQLITE database file which the metadata will be writte into
    args: argparse.Namespace
        Argument used for generating the OPERA database

    '''

    str_datetime_now = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S.%f')
    if not os.path.exists(sql_path):
        raise FileNotFoundError(f'Cannot find SQLITE file: {sql_path}')
    conn = sqlite3.connect(sql_path)
    curs = conn.cursor()

    dict_field_datatype = OrderedDict()
    dict_field_datatype['version'] = ['TEXT', burst_database_version]
    dict_field_datatype['src_database'] = ['TEXT', os.path.basename(args.sqlite_path_in)]
    dict_field_datatype['margin_x'] = ['REAL', args.mxy[0]]
    dict_field_datatype['margin_y'] = ['REAL', args.mxy[1]]
    dict_field_datatype['snap_x'] = ['INTEGER', args.sxy[0]]
    dict_field_datatype['snap_y'] = ['INTEGER', args.sxy[1]]
    dict_field_datatype['datetime_last_mod'] = ['TEXT', str_datetime_now]

    # create the table for metadata inside database
    query_create_table = 'CREATE TABLE IF NOT EXISTS metadata ('
    for field_datatype in dict_field_datatype.items():
        query_create_table += f'{field_datatype[0]} {field_datatype[1][0]}, '
    query_create_table = query_create_table[:-2] + ');'

    curs.execute(query_create_table)

    # Check if there are any records in the metadata table;
    # Clean the table if there any
    curs.execute('SELECT * FROM metadata')
    rec_out = curs.fetchall()
    if len(rec_out) > 0:
        curs.execute('DELETE FROM metadata')

    query_insert = ( 'INSERT INTO metadata '
                    f'{tuple(dict_field_datatype.keys())} '
                     'VALUES('.replace('\'',''))
    for field_item in dict_field_datatype.items():
        if field_item[1][0] == 'TEXT':
            query_insert += f'\"{field_item[1][1]}\", '
        else:
            query_insert += f'{field_item[1][1]}, '
    query_insert = query_insert[:-2] + ')'

    curs.execute(query_insert)
    conn.commit()

    conn.close()


def get_metadata(sql_path:str):
    '''
    Get the metadata (setting, version, creation date, etc.) from db file as dict

    Parameters:
    -----------
    sql_path: str
        Path to the SQLITE database file which the metadata will be writte into
    format: argparse.Namespace

    Return:
    -------
    dict_out: OrderedDict
        metadata

    '''

    if not os.path.exists(sql_path):
        raise FileNotFoundError(f'Cannot find SQLITE file: {sql_path}')
    conn = sqlite3.connect(sql_path)
    curs = conn.cursor()

    curs.execute('SELECT * FROM metadata')
    list_field = [field[0] for field in curs.description]
    records_out = curs.fetchone()


    dict_out = OrderedDict()
    for i_field, field in enumerate(list_field):
        dict_out[field] = records_out[i_field]

    return dict_out
