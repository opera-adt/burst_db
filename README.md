# Burst_DB

Sentinel-1 Burst and Frame Databases for OPERA CSLC-S1/DISP-S1

## Installation

Follow the steps below to install `burst_db` using conda environment.

1. Download source code:

```bash
git clone https://github.com/opera-adt/burst_db
cd burst_db
```

2. Install dependencies:

```bash
conda env create
```

3. Install via pip:

```bash
# run "pip install -e" to install in development mode
python -m pip install .
```

## Usage

Installing the package creates the `opera-db` command line tool:

```bash
$ opera-db --help
Usage: opera-db [OPTIONS] COMMAND [ARGS]...

  Create/interact with OPERA's Sentinel-1 burst/frame databases.

Options:
  --help  Show this message and exit.

Commands:
  create                Generate the OPERA frame database for Sentinel-1...
  historical            Sub-commands for interacting with the historical...
  intersect             Query for frames intersecting a given bounding...
  lookup                Query the geopackage database for one frame ID.
  make-burst-catalog    Create a burst catalog and consistent burst JSON.
  make-reference-dates  Generate a reference dates JSON file for InSAR...
  urls-for-frame        Retrieve URLs for a specific FRAME_ID using the...

```

## Creating the Burst database sqlite file

The `opera-db create` CLI will create the Sqlite database containing the burst IDs, bounding boxes, and UTM EPSG codes for all Sentinel-1 burst ID footprints.

The program uses the database of Sentinel-1 bursts released by ESA. The data can be downloaded from [here](https://sar-mpc.eu/files/S1_burstid_20220530.zip), but if it is not present in the current directory, the program will download it automatically.

A larger GeoPackage is created which contains the burst footprint geometries, which can be viewed/queried with GIS program.


### Frame database information

The other files created from `opera-db create` provide information for the Displacement frame. There are JSON files which map the burst IDs to frame IDs, and frame IDs to burst IDs.

The format of the frame-to-burst mapping is
```python
{
    "data" : {
        "1": {
            "epsg": 32631,
            "is_land": False,
            "is_north_america": False,
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
```python
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

The command also makes a full [Geopackage database](https://www.geopackage.org/) (which is based on sqlite), where the `burst_id_map` table contains the burst geometries, the `frames` table contains the frame geometries, and the `frames_bursts` table is the JOIN table for the many-to-many relationship.
An example SQL query to view all columns of these tables is
```sql
SELECT *
FROM frames f
JOIN frames_bursts fb ON fb.frame_fid = f.fid
JOIN burst_id_map b ON fb.burst_ogc_fid = b.ogc_fid
LIMIT 1;
```
You can also drag the `opera-s1-disp.gpkg` file into QGIS to load the `frames` and `burst_id_map` tables to filter/view the geometries.

## Creating the C
