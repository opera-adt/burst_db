import fnmatch
import zipfile
from pathlib import Path

import geopandas as gpd
import pandas as pd
import requests
import unzip_http
from shapely.geometry import MultiPolygon

USGS_LAND_URL = "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip"  # noqa
GREENLAND_URL = "https://geodata.lib.utexas.edu/download/file/stanford-sd368wz2435-shapefile.zip"  # noqa


def get_usgs_land(outpath=None):
    f"""Get the USGS land data from {USGS_LAND_URL}."""
    outpath = Path(outpath) if outpath else Path.cwd()
    rzf = unzip_http.RemoteZipFile(USGS_LAND_URL)
    # Level 1: Continental land masses and ocean islands, except Antarctica.
    # Level 6: Antarctica based on grounding line boundary.
    paths = ["GSHHS_shp/h/GSHHS_h_L1.*", "GSHHS_shp/h/GSHHS_h_L6.*"]
    shp_files = []
    dfs = []
    for fn in rzf.infolist():
        if not any(fnmatch.fnmatch(fn.filename, g) for g in paths):
            continue
        outname = outpath / fn.filename
        if outname.suffix == (".shp"):
            shp_files.append(outname)
        if not outname.exists():
            outname.parent.mkdir(parents=True, exist_ok=True)
            with rzf.open(fn) as fp, open(outname, "wb") as fout:
                print(f"Extracting {fn.filename} to {outname}")
                while r := fp.read(2**18):
                    fout.write(r)
    for p in shp_files:
        dfs.append(gpd.read_file(p))
    return dfs


def get_land_df(
    buffer_deg=0.2,
    outname="usgs_land_{d}deg_buffered.geojson",
    driver="GeoJSON",
    zip=True,
):
    outname = outname.format(d=buffer_deg)
    if outname and Path(outname).exists():
        print(f"Loading {outname} from disk")
        return gpd.read_file(outname)
    elif Path(outname + ".zip").exists():
        print(f"Loading {outname}.zip from disk")
        return gpd.read_file(str(outname) + ".zip")

    # If we haven't already made the file, make it
    df_land_cont, df_antarctica = get_usgs_land()
    df_land = pd.concat([df_land_cont, df_antarctica], axis=0)[["geometry"]]
    df_land.geometry = df_land.geometry.buffer(buffer_deg)
    df_land = df_land.dissolve()

    df_land.to_file(outname, driver=driver)
    if zip and outname.endswith(".geojson"):
        outname_zipped = Path(str(outname) + ".zip")
        # zip and remove the original
        with zipfile.ZipFile(
            outname_zipped, "w", compression=zipfile.ZIP_DEFLATED
        ) as zf:
            zf.write(outname)

        # Remove the original
        Path(outname).unlink()

    return df_land


def get_greenland_shape(outpath=None, buffer_deg=0.2) -> MultiPolygon:
    f"""Get the Greenland data from {GREENLAND_URL}."""
    outpath = Path(outpath) if outpath else Path.cwd()
    outname = outpath / f"greenland_{buffer_deg}deg_buffered.geojson"
    if outname.exists():
        print(f"Loading {outname} from disk")
        return gpd.read_file(outname).iloc[0].geometry

    # download the whole greenland shapefile
    print("Downloading Greenland shapefile...")
    r = requests.get(GREENLAND_URL)
    zipfile = outpath / "greenland.zip"
    if not zipfile.exists():
        with open(zipfile, "wb") as fout:
            fout.write(r.content)

    df = gpd.read_file(zipfile)
    print("Simplifying and buffering Greenland shapefile...")
    g = df.iloc[0].geometry
    # Now simplify first, then buffer
    gs = g.simplify(0.1)
    g_buffered = gs.buffer(buffer_deg)
    # Save for later
    gpd.GeoDataFrame(geometry=[g_buffered]).to_file(outname, driver="GeoJSON")
