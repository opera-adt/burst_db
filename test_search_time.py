'''
A script to test the search performance of each data format and loading scmheme
(CSV, JSONized CSV, JSON, SQLITE)

'''

import json
import os
import random
import time

import sqlite3


def get_records(list_burst_id:list, path_table:str):
    '''Placeholder'''
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

    elif path_table.lower().endswith('.json'):
        with open(path_table, 'r', encoding='utf-8') as fin_json:
            dict_burst = json.load(fin_json)

        # Perform searching for the burst information
        for burst_id in list_burst_id:
            _ = dict_burst[burst_id]


    elif path_table.lower().endswith('.sqlite3'):
        with sqlite3.connect(path_table) as conn_select:
            cur_select = conn_select.cursor()
            for burst_id in list_burst_id:
                str_sql_select = f'''SELECT * FROM burst WHERE burst_id="{burst_id}"'''
                cur_select.execute(str_sql_select)
                _ = cur_select.fetchall()
    else:
        raise ValueError(f'Data file not supported: {path_table}')



def get_records_old_csv(list_burst_id:str, path_table:str):
    '''Find the burst coverage information from the CSV table'''
    with open(path_table, 'r', encoding='utf-8') as fin_csv:
        lines_csv = fin_csv.readlines()
        if lines_csv[-1] == '':  # Trim the last white line when necessary
            lines_csv = lines_csv[:-1]

    # Convert CSV to dict
    num_record = len(lines_csv) - 1
    list_burst_id = [None]*num_record

    for i_record in range(1, num_record+1):
        token_line = lines_csv[i_record].split(', ')
        list_burst_id[i_record-1] = token_line[0]

    # Perform searching for the burst information
    for burst_id in list_burst_id:
        _ = list_burst_id.index(burst_id)


if __name__=='__main__':
    PATH_REPO = os.path.dirname(__file__)
    CSV_OUT = f'{PATH_REPO}/output/Burst_Coverage_Metadata.csv'
    JSON_OUT = f'{PATH_REPO}/output/Burst_Coverage_Metadata.json'
    SQLITE_OUT = f'{PATH_REPO}/output/Burst_Coverage_Metadata.sqlite3'

    # load the list of burst id
    with open(CSV_OUT, 'r', encoding='utf-8') as fin:
        lines_csv_in = fin.readlines()
        if lines_csv_in[-1] == '':  # Trim the last white line when necessary
            lines_csv_in = lines_csv_in[:-1]

    num_record_csv = len(lines_csv_in)-1
    list_burst_id_query = [None] * num_record_csv
    for i in range(1, num_record_csv+1):
        token_line_csv = lines_csv_in[i].split(', ')
        list_burst_id_query[i-1] = token_line_csv[0]

    print('Populating burst IDs to search')
    NUM_SEARCH = 1000000
    list_burst_id_search = [None] * NUM_SEARCH
    for i in range(NUM_SEARCH):
        id_put = random.randint(0,num_record_csv-1)
        list_burst_id_search[i] = list_burst_id_query[id_put]


    print('Testing CSV...', end='')
    t0_csv = time.time()
    get_records(list_burst_id_search, CSV_OUT)
    t1_csv = time.time()
    print('COMPLETED.')

    print('Testing CSV Old fashion...', end='')
    t0_csv_old = time.time()
    get_records_old_csv(list_burst_id_search, CSV_OUT)
    t1_csv_old = time.time()
    print('COMPLETED.')

    print('Testing JSON...', end='')
    t0_json = time.time()
    get_records(list_burst_id_search, JSON_OUT)
    t1_json = time.time()
    print('COMPLETED.')

    print('Testing SQLite...', end='')
    t0_sqlite = time.time()
    get_records(list_burst_id_search, SQLITE_OUT)
    t1_sqlite = time.time()
    print('COMPLETED.')


    print('RESULTS:')
    print('-------------------------------------------')
    print(f'CSV    : {t1_csv - t0_csv : 0.6f} seconds')
    print(f'JSON   : {t1_json - t0_json : 0.6f} seconds')
    print(f'SQLite : {t1_sqlite - t0_sqlite : 0.6f} seconds')
    print(f'CSV OLD: {t1_csv_old - t0_csv_old : 0.6f} seconds')
