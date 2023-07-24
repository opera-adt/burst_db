# Burst_DB
Sentinel-1 Burst Coverage Database for OPERA SAS

🚨 This toolbox is still in pre-alpha stage and undergoing rapid development. 🚨

## How to install
Follow the steps below to install `burst_db` using conda environment.

1. Download source code:

```bash
git clone https://github.com/opera-adt/burst_db burst_db
cd burst_db
```

2. Install dependencies:

```bash
conda install -c conda-forge --file environment.yml
```

3. Install `burst_db` via pip:

```bash
# run "pip install -e" to install in development mode
python -m pip install .
```

## How to use
`build_database.py` provides CLI to provide the input / output files and settings for the bounding box calculation. Below is the usage of the program.

>build_database.py [-h] [-mxy MXY MXY] [-sxy SXY SXY] [-d DEPLOYABLE] sqlite_path_in sqlite_path_out

- `-mxy` : x / y margin to apply to the bounding box. Default: [5000.0 5000.0]
- `-sxy` : x / y snap value to ceil/floor the boundaries' coordinates. Default: [30 30]
- `-d` : (Optional) Smaller version of the DB will be saved to the path `DEPLOYABLE`
- `sqlite_path_in` : Path to the source SQLite database file for Sentinel-1 burst map released by ESA. The data can be downloaded from [here](https://sar-mpc.eu/files/S1_burstid_20220530.zip).
- `sqlite_path_out` : Path to the output SQLite database file.


## Frame database information

After running `pip install .` , the `opera-create-db` command will create the sqlite Frame Database, as well as JSON files which map the burst IDs to frame IDs, and frame IDs to burst IDs.

The format of the frame-to-burst mapping is
```json
{
    "data" : {
        "1": {
            "epsg": 32631,
            "is_land": False,
            "xmin": 500160,
            "ymin": 78240,
            "xmax": 789960,
            "ymax": 322740,
            "burst_ids": [
                "t001_000001_iw1",
                "t001_000001_iw2",
                "t001_000001_iw3",
                "t001_000002_iw1",
                ...
                "t001_000009_iw3"
            ]
        }, ...
    },
    "metadata": {
        "version": "0.1.2", "margin": 5000.0, ...
    }
}
```
where the keys of the the `data` dict are the frame IDs.

The burst-to-frame mapping has the structure
```json
{
    "data" : {
        "t001_000001_iw1": {"frame_ids": [1]},
        "t001_000001_iw2": {"frame_ids": [1]},
        ...
    },
    "metadata": {
        "version": "0.1.2", "margin": 5000.0, ...
    }
}
```
These data structures can be read into python using the function `build_frame_db.read_zipped_json` .

The `opera-create-db` command also makes the full [Geopackage database](https://www.geopackage.org/) (which is based on sqlite), where the `burst_id_map` table contains the burst geometries, the `frames` table contains the frame geometries, and the `frames_bursts` table is the JOIN table for the many-to-many relationship.
An example SQL query to view all columns of these tables is
```sql
SELECT *
FROM frames f
JOIN frames_bursts fb ON fb.frame_fid = f.fid
JOIN burst_id_map b ON fb.burst_ogc_fid = b.ogc_fid
LIMIT 1;
```
You can also drag the `opera-s1-disp.gpkg` file into QGIS to load the `frames` and `burst_id_map` tables to filter/view the geometries.


### License
**Copyright (c) 2022** California Institute of Technology (“Caltech”). U.S. Government
sponsorship acknowledged.

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided
that the following conditions are met:
* Redistributions of source code must retain the above copyright notice, this list of conditions and
the following disclaimer.
* Redistributions in binary form must reproduce the above copyright notice, this list of conditions
and the following disclaimer in the documentation and/or other materials provided with the
distribution.
* Neither the name of Caltech nor its operating division, the Jet Propulsion Laboratory, nor the
names of its contributors may be used to endorse or promote products derived from this software
without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
