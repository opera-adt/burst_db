# Burst_DB
Sentinel-1 Burst Coverage Database for OPERA SAS

🚨 This toolbox is still in pre-alpha stage and undergoing rapid development. 🚨

## How to install
Follow the steps below to install `burst_db` using conda environment.

1. Download source code:

```bash
git clone https://github.com/opera-adt/burst_db burst_db
```

2. Install dependencies:

```bash
conda install -c conda-forge --file burst_db/requirements.txt
```

3. Install `burst_db` via pip:

```bash
# run "pip install -e" to install in development mode
python -m pip install ./burst_db
```

## How to use
`build_database.py` provides CLI to provide the input / output files and settings for the bounding box calculation. Below is the usage of the program.

>build_database.py [-h] [-mxy MXY MXY] [-sxy SXY SXY] [-d DEPLOYABLE] sqlite_path_in sqlite_path_out

- `-mxy` : x / y margin to apply to the bounding box. Default: [5000.0 5000.0]
- `-sxy` : x / y snap value to ceil/floor the boundaries' coordinates. Default: [30 30]
- `-d` : (Optional) Smaller version of the DB will be saved to the path `DEPLOYABLE`

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
