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

## Creating a new release

After making changes to the code, a new release can be created by running the following commands:

```bash
# For example, if the new version is 0.12.0
git tag v0.12.0
pip install -e .

# Setup in a new folder
mkdir -p test_012 && cd test_012
# Copy last CMR survey of CSLC products
cp ../test_011/cmr_survey.2016-07-01_to_2024-12-31.2025-06-12.opera-pcm-3.1.6.csv.tar.gz .
# OR: if updating, get the new survey from SDS
make -f ../Makefile
```

The result will be a folder with the following files:

```bash
$ ls test_012
burst_map_IW_000001_375887.sqlite3                                          opera-disp-s1-blackout-dates-2025-08-12.json
burst-id-geometries-simple-0.12.0.geojson                                   opera-disp-s1-consistent-burst-ids-2025-08-12-2016-07-01_to_2024-12-31.json
burst-id-geometries-simple-0.12.0.geojson.zip                               opera-disp-s1-consistent-burst-ids-no-blackout.json
cmr_survey_2016-07-01_to_2024-12-31.csv                                     opera-disp-s1-reference-dates-2025-08-12.json
cmr_survey.2016-07-01_to_2024-12-31.2025-06-12.opera-pcm-3.1.6.csv.tar.gz   opera-s1-disp-0.12.0-2d.gpkg
cslc-burst-database-2025-08-12.duckdb                                       opera-s1-disp-0.12.0-burst-to-frame.json.zip
frame-geometries-simple-0.12.0.geojson                                      opera-s1-disp-0.12.0-frame-to-burst.json.zip
frame-geometries-simple-0.12.0.geojson.zip                                  opera-s1-disp-0.12.0.gpkg
GSHHS_shp                                                                   usgs_land_0.3deg_buffered.geojson.zip
opera-burst-bbox-only.sqlite3
```

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

## Creating the databases

![production flow](./docs/DISP-S1-database-production-flow.drawio.svg)

The following instructions show how to create the auxilliary database for DISP-S1 processing.

The example uses version 0.7.0.

### Prerequisites

- CMR survey file (gzipped)
- Snow blackout dates file (optional)

### Steps

1. Set up the output folder for the current release:

```bash
mkdir outputs-070
cd outputs-070
```

2. Copy in the CMR survey gzipped file and (optionally) snow blackout file

```bash
cp /path/to/survey/cmr-surveys/cmr_survey.2016-07-01_to_2024-12-10.csv.tar.gz .
cp /path/to/opera-disp-s1-blackout-dates-2024-10-16.json .
```

3. Run `make`

```bash
make -f ../Makefile
```

Typical processing should take ~5-8 minutes, depending on download speed.

### Note on optional prerequisite input: Snow blackout dates

To create the blackout input, the following module is used (currently no CLI):

```python
import burst_db.create_blackout_dates_s1
burst_db.create_blackout_dates_s1.gdf_to_blackout_json(input_json_file)
```

### Details on subcommands run

#### Frame ID and Burst ID Geopackage database

```bash
opera-db create # Creates the geopackage, and aux. geojson helper files
```

Parse the "CMR survey" of all existing bursts, and keep a set which is consistent through space and time (i.e. no spatial gaps will appear while processing)

```bash
opera-db make-burst-catalog ...
```

Set up a JSON, one key per DISP-S1 Frame ID, listing the "reference date changes".
This indicates to the processing system that we should start outputting data with respect to a new reference, to avoid attempting to form very long temporal baseline interferograms.

```bash
opera-db make-reference-dates
```

## Reconciling and Labeling Burst Databases

The `reconcile_and_label_db.py` script is a standalone tool for reconciling differences between two burst database JSON files and adding processing mode labels to sensing times.

### Features

1. **Database Reconciliation**: Compares old and new consistent burst database JSON files and reconciles differences:
   - If the new database has more burst IDs than the old for a frame, it uses the old burst IDs
   - If the new database is missing sensing times from the old, it adds them back
   - If there's no overlap in sensing times (indicating a complete restart after a gap), the new data is kept as-is

2. **Processing Mode Labeling**: Assigns processing mode labels to each sensing time:
   - `historical_XX`: Full batches of 15 sensing times (ready for historical processing)
   - `forward_XX`: Partial batches with fewer than 15 sensing times (for forward processing)
   - `no_run`: Groups with fewer than 15 total sensing times (insufficient data)
   - The `_XX` suffix (e.g., `_01`, `_02`) indicates the group number, which increments after temporal gaps of 2+ years

### Usage

```bash
# Reconcile two databases and add processing mode labels
python src/burst_db/reconcile_and_label_db.py \
    --old-db old_consistent_burst_ids.json \
    --new-db new_consistent_burst_ids.json \
    --output labeled_output.json

# Only add processing mode labels (skip reconciliation)
python src/burst_db/reconcile_and_label_db.py \
    --new-db input.json \
    --output output.json \
    --no-reconcile

# Customize batch size and gap threshold
python src/burst_db/reconcile_and_label_db.py \
    --old-db old.json \
    --new-db new.json \
    --output output.json \
    --batch-size 20 \
    --gap-threshold 1.5 \
    --verbose
```

### Command Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--old-db` | Path to the old burst database JSON file | Required (unless `--no-reconcile`) |
| `--new-db` | Path to the new burst database JSON file | Required |
| `--output` | Path for the output JSON file | Required |
| `--no-reconcile` | Skip reconciliation, only add labels | False |
| `--batch-size` | Number of sensing times per batch | 15 |
| `--gap-threshold` | Gap threshold in years to restart batching | 2.0 |
| `--verbose` | Print detailed frame information | False |

### Output Format

The output JSON replaces the `sensing_time_list` array with a dictionary mapping sensing times to their labels:

```json
{
  "metadata": { ... },
  "data": {
    "831": {
      "burst_id_list": ["t004_006645_iw1", "t004_006646_iw1", ...],
      "sensing_time_list": {
        "2016-07-02T23:05:35": "historical_01",
        "2016-09-24T23:05:39": "historical_01",
        ...
        "2025-10-19T23:06:08": "forward_01"
      }
    },
    ...
  }
}
```
