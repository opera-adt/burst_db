#!/usr/bin/env python
'''
A script to test the search performance of each data format and loading scmheme
(CSV, JSONized CSV, JSON, SQLITE)

'''

import argparse
import csv
import json
import random
import time

import sqlite3
import numpy as np



def populate_bursts_to_search(path_data:str, num_query:int):
    '''
    Randomly select the burst IDs. Return the selections as list
    Supported file format: CSV, JSON, sqlite3.

    Parameters:
    -----------
    path_data: str
        Path to the input burst database (CSV, JSON, SQLITE, SQLITE3)

    num_query: int
        Length of the output list of the burst

    Return:
    -------
    list_burst_id_search: list
        List of the burst ID to query
    '''

    if path_data.lower().endswith('.sqlite') or path_data.lower().endswith('.sqlite3'):
        conn = sqlite3.connect(path_data)
        curs = conn.cursor()

        query_out = curs.execute('SELECT burst_id from burst')
        records_out = query_out.fetchall()
        list_burst_all = [record[0] for record in records_out]

        curs.close()
        conn.close()

    elif path_data.lower().endswith('.csv'):
        with open(path_data, 'r', encoding='utf-8') as fin:
            reader_csv = csv.reader(fin)
            list_burst_all = [row[0] for row in reader_csv]

    elif path_data.lower().endswith('.json'):
        with open(path_data, 'r', encoding='utf-8') as fin:
            json_in = json.load(fin)
        list_burst_all = list(json_in.keys())

    else:
        raise ValueError(f'Cannot extract the bursts to search from: {path_data}')

    # Randonly extract the bursts to find
    num_burst_all = len(list_burst_all)
    list_burst_id_search = [None] * num_query
    for id_search in range(num_query):
        id_put = random.randint(0,num_burst_all-1)
        list_burst_id_search[id_search] = list_burst_all[id_put]

    return list_burst_id_search


def query_burst(list_burst_to_find: list, path_data: str, old_csv=False):
    '''
    Determine what query scheme to use

    Parameters:
        list_burst_to_find: list
            List of burst IDs to search from `path_data`
        path_data: str
            Path to the burst database
        old_csv: bool
            If True, CSV will be loaded and treated as a dict.
            If False, CSV will be loaded and treated as a list

    '''
    if path_data.lower().endswith('.sqlite') or path_data.upper().endswith('.sqlite3'):
        query_burst_sqlite(list_burst_to_find, path_data)
    elif path_data.lower().endswith('.csv'):
        if old_csv:
            query_burst_old_csv(list_burst_to_find, path_data)
        else:
            query_burst_csv(list_burst_to_find, path_data)
    elif path_data.lower().endswith('.json'):
        query_burst_json(list_burst_to_find, path_data)
    else:
        raise ValueError(f'Cannot determine the format of data file: {path_data}')


def query_burst_csv(list_burst_id:list, path_table:str):
    '''
    Find the records from CSV using burst_id.
    This function internally converts CSV into dict, and search the records inside dict.

    Parameters:
    -----------
    list_burst_id: list
        List of burst IDs to query

    path_table: str
        Path to the .csv file

    '''
    if path_table.lower().endswith('.csv'):
        with open(path_table, 'r', encoding='utf-8') as fin_csv:
            lines_csv = fin_csv.readlines()
            if lines_csv[-1] == '':  # Trim the last white line when necessary
                lines_csv = lines_csv[:-1]

        # Convert CSV to dict
        num_record = len(lines_csv)-1
        dict_burst = {}
        for i_record in range(1,num_record+1):
            token_line = lines_csv[i_record].split(', ')
            dict_burst[token_line[0]]={
                'EPSG':int(token_line[1]),
                'xmin':float(token_line[2]),
                'ymin':float(token_line[3]),
                'xmax':float(token_line[4]),
                'ymax':float(token_line[5])
            }

        for burst_id in list_burst_id:
            _ = dict_burst[burst_id]

    else:
        raise ValueError(f'Data file not supported: {path_table}')


def query_burst_json(list_burst_id:list, path_table:str):
    '''
    Find the records from JSON file using burst_id.

    Parameters:
    -----------
    list_burst_id: list
        List of burst IDs to query

    path_table: str
        Path to the .json file

    '''

    if path_table.lower().endswith('.json'):
        with open(path_table, 'r', encoding='utf-8') as fin_json:
            dict_burst = json.load(fin_json)

        # Perform searching for the burst information
        for burst_id in list_burst_id:
            _ = dict_burst[burst_id]

    else:
        raise ValueError(f'Data file not supported: {path_table}')


def query_burst_sqlite(list_burst_id:list, path_table:str):
    '''
    Find the records from sqlite database using burst_id.

    Parameters:
    -----------
    list_burst_id: list
        List of burst IDs to query

    path_table: str
        Path to the .sqlite or .sqlite3 file

    '''

    if path_table.lower().endswith('.sqlite'):
        with sqlite3.connect(path_table) as conn_select:
            cur_select = conn_select.cursor()
            for burst_id in list_burst_id:
                str_sql_select = f'''SELECT * FROM burst WHERE burst_id="{burst_id}"'''
                cur_select.execute(str_sql_select)
                _ = cur_select.fetchall()
    else:
        raise ValueError(f'Data file not supported: {path_table}')


def query_burst_old_csv(list_burst_id_search:list, path_table:str):
    '''
    Find the records from CSV using burst_id.
    This function internally converts CSV into list, and search the records inside the list.

    Parameters:
    -----------
    list_burst_id: list
        List of burst IDs to query

    path_table: str
        Path to the .csv file

    '''

    with open(path_table, 'r', encoding='utf-8') as fin_csv:
        lines_csv = fin_csv.readlines()

        # Trim the last white line when necessary
        if lines_csv[-1] == '':
            lines_csv = lines_csv[:-1]

    # Convert CSV to dict
    num_record = len(lines_csv) - 1
    list_burst_id = [None]*num_record

    for i_record in range(1, num_record+1):
        token_line = lines_csv[i_record].split(', ')
        list_burst_id[i_record-1] = token_line[0]

    # Perform searching for the burst information
    for burst_id_search in list_burst_id_search:
        _ = list_burst_id.index(burst_id_search)


if __name__=='__main__':
    parser = argparse.ArgumentParser(
        description='Performance tester for Sentinel-1 Burst coveage database'
    )

    parser.add_argument('-n','--num_query',
                        type=int, default=30, help='Number of queries per round. Default = 30')
    parser.add_argument('-r','--repeat',
                        type=int, default=20, help='Number of rounds. Default = 20')
    parser.add_argument('--old_csv',
                        default=False, action="store_true",
                        help='option for .csv data to treat the dta file as list. Default = False')
    parser.add_argument('path_data',
                        type=str, help='Path to the data file (sqlite, sqlite3, csv, json')

    args=parser.parse_args()


    # printout the user input
    print('\nTest setting:')
    print('-------------')
    print(f'path_data                  : {args.path_data}')
    print(f'Number of bursts per round : {args.num_query}')
    print(f'Number of rounds           : {args.repeat}')
    print(f'old_csv                    : {args.old_csv}\n')


    # Numpy array to keep record of the processing time for each round
    list_time = np.zeros(args.repeat)
    for i_repeat in range(args.repeat):
        print(f'Round {i_repeat+1} / {args.repeat} -',end=' ')
        # Populate the data to search
        list_burst_query = populate_bursts_to_search(args.path_data, args.num_query)

        t0 = time.time()
        query_burst(list_burst_query, args.path_data, args.old_csv)
        t1 = time.time()

        list_time[i_repeat] = t1 - t0
        print(f'{list_time[i_repeat]} sec.        ',end='\r')

    print('\nTest complete:')
    print(f'min   : {list_time.min():0.6f} sec. '
          f'({list_time.min()/args.num_query*1000.0:0.6f} ms. per burst)')
    print(f'max   : {list_time.max():0.6f} sec. '
          f'({list_time.max()/args.num_query*1000.0:0.6f} ms. per burst)')
    print(f'mean  : {list_time.mean():0.6f} sec. '
          f'({list_time.mean()/args.num_query*1000.0:0.6f} ms. per burst)')
    print(f'stdev : {list_time.std():0.6f} sec. '
          f'({list_time.std()/args.num_query*1000.0:0.6f} ms. per burst)')
