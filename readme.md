# Sentinel-1 Burst Coverage Database


## Usage

To create a database of the Sentinel-1 burst footprints, use the `make_db.py` script. 
This script will download ESA's burst map (available here: https://sar-mpc.eu/files/S1_burstid_20220530.zip) and augment it with a unique ID per burst, as well as the bounding box limits in UTM. 


```console
$ python make_db.py --output-path burst_map.sqlite3

ESA database not found, downloading...
Downloading ESA burst database
...
```

If you already have the ESA database, you can use the `--esa-db` option to specify the path to the database. 
Otherwise, `wget` is used in the script to download the files.


The full set of options are available through the `--help` flag:

```console
$ python make_db.py --help
usage: make_db.py [-h] [--esa-db-path ESA_DB_PATH] [--output-path OUTPUT_PATH] [--snap SNAP] [--margin MARGIN] [--limit LIMIT] [--max-procs MAX_PROCS]

options:
  -h, --help            show this help message and exit
  --esa-db-path ESA_DB_PATH
                        Path to the ESA sqlite burst database to convert, downloaded from https://sar-mpc.eu/files/S1_burstid_20220530.zip . Will be downloaded if not exists. (default: burst_map_IW_000001_375887.sqlite3)
  --output-path OUTPUT_PATH
                        Path to the output database (default: burst_map_IW_000001_375887.OPERA-JPL.20221109_200032.sqlite3)
  --snap SNAP           Snap the bounding box to the nearest multiple of this value. (default: 50.0)
  --margin MARGIN       Add this margin surrounding the bounding box of bursts. (default: 4000.0)
  --limit LIMIT         For testing, limit the number of rows to process (default: None)
  --max-procs MAX_PROCS
                        Max CPU count to use for processing geometries (default: None)
```

## Setup

```console
# Download the code
git clone git@github.com:opera-adt/burst_db.git
# Install the requirements
conda create -n my-burst-env --file requirements.txt
conda activate my-burst-env

python make_db.py --help


```

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
