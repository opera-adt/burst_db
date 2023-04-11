import fnmatch
from pathlib import Path

import geopandas as gpd
import pandas as pd
import unzip_http

USGS_LAND_URL = (
    "https://www.ngdc.noaa.gov/mgg/shorelines/data/gshhg/latest/gshhg-shp-2.3.7.zip"
)


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
            with rzf.open(fn) as fp, open(outname, "wb") as fout:
                print(f"Extracting {fn.filename} to {outname}")
                while r := fp.read(2**18):
                    fout.write(r)
    for p in shp_files:
        dfs.append(gpd.read_file(p))
    return dfs


def get_land_df(
    buffer_deg=0.2, outname="usgs_land_{d}deg_buffered.geojson", driver="GeoJSON"
):
    outname = outname.format(d=buffer_deg)
    if outname and Path(outname).exists():
        print(f"Loading {outname} from disk")
        return gpd.read_file(outname)

    df_land_cont, df_antarctica = get_usgs_land()
    df_land = pd.concat([df_land_cont, df_antarctica], axis=0)[["geometry"]]
    df_land.geometry = df_land.geometry.buffer(buffer_deg)
    df_land = df_land.dissolve()
    if outname:
        df_land.to_file(outname, driver=driver)
    return df_land
