import sqlite3
import json
import os

path_home = os.getenv('HOME')

PATH_AUGMENTED_BURST_MAP = (f'{path_home}/Documents/DATA/Sensor/Sentinel-1/'
                            'S1_burstid_20220530/IW/sqlite/'
                            'burst_map_IW_000001_375887.OPERA-JPL.20220818_113416.sqlite3')

PATH_TOP=os.path.dirname(__file__)
PATH_CSV_OUT = f'{PATH_TOP}/data/CSV_OUT_ALL_BURSTS.csv'
PATH_JSON_OUT = f'{PATH_TOP}/data/CSV_OUT_ALL_BURSTS.json'
PATH_SQLITE_OUT = f'{PATH_TOP}/data/CSV_OUT_ALL_BURSTS.sqlite'

conn=sqlite3.connect(PATH_AUGMENTED_BURST_MAP)
curs=conn.cursor()
curs.execute('SELECT burst_id_jpl, EPSG, xmin, ymin, xmax, ymax FROM burst_id_map')
query_out=curs.fetchall()
num_record=len(query_out)
curs.close()
conn.close()

#export to CSV
print('Exporting to CSV')
with open(PATH_CSV_OUT,'w+',encoding='utf-8') as fout:
    #write CSV header
    fout.write('burst_id, EPSG, xmim, ymin, xmax, ymax\n')

    for i, record in enumerate(query_out):
        fout.write(f'{record[0]}, {record[1]}, {record[2]}, {record[3]}, {record[4]}, {record[5]}\n' )

# Export to json
print('Exporting to JSON')
dict_export={}
for i in range(num_record):
    dict_export[query_out[i][0]]={
        'EPSG':query_out[i][1],
        'extent':[query_out[i][2],query_out[i][3],query_out[i][4],query_out[i][5]]
    }

with open(PATH_JSON_OUT, 'w+', encoding='utf8') as fout:
    json.dump(dict_export, fout, indent=2)


# Export to SQLITE
print('Exporting to SQLITE')
with sqlite3.connect(PATH_SQLITE_OUT) as conn_out:
    curs_out=conn_out.cursor()
    curs_out.execute('''CREATE TABLE IF NOT EXISTS burst (
        burst_id text PRIMARY KEY,
        EPSG integer,
        xmin float,
        ymin float,
        xmax float,
        ymax float);
        ''')

    for i_burst,record in enumerate(query_out):
        print(f'{i_burst:,} / {num_record:,}',end='\r')
        str_sql_command=f'INSERT INTO burst (burst_id ,EPSG, xmin, ymin, xmax, ymax) '\
                        f'VALUES("{record[0]}", {record[1]}, '\
                        f'{record[2]}, {record[3]}, '\
                        f'{record[4]}, {record[5]})'

        curs_out.execute(str_sql_command)
    print('\n')
    curs_out.execute('CREATE INDEX index_burst ON burst (burst_id)')
    conn_out.commit()
