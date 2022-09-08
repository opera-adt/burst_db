#!/usr/bin/env python
'''
A script to export the burst geogrid from the augmented burst map

'''
import argparse
import sqlite3
import json
import os


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
    curs.execute('SELECT burst_id_jpl, EPSG, xmin, ymin, xmax, ymax FROM burst_id_map')
    records_out = curs.fetchall()

    curs.close()
    conn.close()

    return records_out


def export_to_csv(records_out: list, path_csv_out: str):
    '''
    Write out the burst geogrid information as .csv file

    Parameters:
        records_out: list
            List of burst geogrid information retrieved from `extract_burst_geogrid_data()`
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
            List of burst geogrid information retrieved from `extract_burst_geogrid_data()`
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


def export_to_sqlite(records_out: list, path_sqlite_out: str):
    '''
    Write out the burst geogrid information as sqlite database file

    Parameters:
        records_out: list
            List of burst geogrid information retrieved from `extract_burst_geogrid_data()`
        path_sqlite_out: str
            path to the .sqlite file to write out

    '''

    num_record = len(records_out)

    # Export to SQLITE
    print('Exporting to SQLITE')
    with sqlite3.connect(path_sqlite_out) as conn_out:
        curs_out = conn_out.cursor()
        curs_out.execute('CREATE TABLE IF NOT EXISTS burst ('
                         'burst_id text PRIMARY KEY, EPSG integer, '
                         'xmin float, ymin float, xmax float, ymax float);')

        for i_record, record in enumerate(records_out):
            print(f' Processing: {i_record:,} / {num_record:,}      ', end='\r')
            str_sql_command=f'INSERT INTO burst (burst_id ,EPSG, xmin, ymin, xmax, ymax) '\
                            f'VALUES("{record[0]}", {record[1]}, '\
                            f'{record[2]}, {record[3]}, '\
                            f'{record[4]}, {record[5]})'

            curs_out.execute(str_sql_command)
        print('\n')
        curs_out.execute('CREATE INDEX index_burst ON burst (burst_id)')
        conn_out.commit()


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Extracts burst geogrid data from augmented burstmap. Writes out the data into .csv, .json, or .sqlite')

    parser.add_argument('path_augmented_burstmap', type=str,
                        help='Path to the augmented burst map from which the burst geogrid will be extracted')
    parser.add_argument('path_databse_out', type=str,
                        help='Path to the output database file\nSupported format: .csv, .json, .sqlite')

    args=parser.parse_args()

    if not os.path.exists(args.path_augmented_burstmap):
        raise ValueError(f'Cannot find the augmented burst map: {args.path_augmented_burstmap}')

    records_burst_data = extract_burst_geogrid_data(args.path_augmented_burstmap)

    if parser.path_database_out.endswith('.csv'):
        export_to_csv(records_burst_data, args.path_database_out)
    elif parser.path_database_out.endswith('.json'):
        export_to_json(records_burst_data, args.path_database_out)
    elif (parser.path_database_out.endswith('.sqlite')) or (args.path_database_out.endswith('.sqlite3')):
        export_to_sqlite(records_burst_data, args.path_database_out)
    else:
        raise ValueError(f'Output file format was not supported: {args.path_database_out}')
