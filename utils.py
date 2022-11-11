from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
from pyproj import Transformer
from shapely import ops, wkb

TRANSFORMERS = {}


def get_epsg(geom):
    if len(geom.geoms) == 1:
        point = geom.centroid
    else:
        # antimeridian case
        # Transform to the first UTM zone, get the centroid, then transform back
        # Doesnt matter if it's the right EPSG; just need it to be in
        # projected coordinates to make the centroid easy
        utm_epsg = 32601 if geom.centroid.y > 0 else 32701
        src, dst = 4326, utm_epsg
        point = transform(transform(geom, src, dst).centroid, dst, src)

    return get_point_epsg(point.y, point.x)


def get_utm_bbox(geom, epsg, snap=50.0, margin=1000.0):
    geom_transformed = transform(geom, 4326, epsg)
    return snap_bbox(geom_transformed.bounds, margin, snap)


def transform(geom, src_epsg, dst_epsg):
    if (src_epsg, dst_epsg) in TRANSFORMERS:
        t = TRANSFORMERS[(src_epsg, dst_epsg)]
    else:
        t = Transformer.from_crs(
            src_epsg,
            dst_epsg,
            always_xy=True,
        )
        TRANSFORMERS[(src_epsg, dst_epsg)] = t
    return ops.transform(t.transform, geom)


def to_2d(shape):
    return ops.transform(lambda *args: args[:2], shape)


def _convert_wkb_to_2d(b):
    return to_2d(wkb.loads(b))
    # return wkb.dumps(to_2d(wkb.loads(b)))


def _process_geom(geom_wkb, snap=50.0, margin=1000.0):
    geom_2d = _convert_wkb_to_2d(geom_wkb)
    epsg = get_epsg(geom_2d)
    bbox = get_utm_bbox(geom_2d, epsg, snap=snap, margin=margin)
    return wkb.dumps(geom_2d), epsg, *bbox


def process_geom_parallel(wkb_series, max_procs=30, snap=50.0, margin=1000.0):
    """Loop through the burst geometries and find the EPSG/bounding boxes.

    Consists of 3 parts:
    1. Convert the WKB to 2D
    2. Find the EPSG code for each burst using the centroid
    3. Transform the bounding box for each burst to UTM coordinates

    Parameters
    ----------
    wkb_series : Iterable[bytes]
        Iterable of WKB geometries
    max_procs : int, optional
        Maximum number of processes to use.
        By default 30, or the number of CPUs, whichever is lower.
    snap : float, optional
        Snap the bounding box to the nearest multiple of this value.
        By default 50.0 meters.
    margin : float, optional
        Add this margin to the bounding box.
        By default 1000.0 meters.


    Returns
    -------
    geoms_2d : List[bytes]
        List of WKB geometries in 2D
    epsg_codes : List[int]
        List of EPSG codes
    bboxes : List[Tuple[float]]
        List of bounding boxes
    """
    workers = min(max_procs, cpu_count())
    # return list(map(_process_geom, wkb_series))
    print(f"Using {workers} cores to process geometries")
    func = partial(_process_geom, snap=snap, margin=margin)
    with Pool(workers) as p:
        # passing through kwargs is annoying
        return list(
            p.map(
                func,
                wkb_series,
            )
        )


def snap_bbox(bbox, margin=1000, snap=50):
    xmin, ymin, xmax, ymax = bbox

    if margin > 0:
        xmin -= margin
        xmax += margin
        ymin -= margin
        ymax += margin

    if snap > 0:
        xmin = int(np.floor(xmin / snap) * snap)
        xmax = int(np.ceil(xmax / snap) * snap)
        ymin = int(np.floor(ymin / snap) * snap)
        ymax = int(np.ceil(ymax / snap) * snap)

    return (xmin, ymin, xmax, ymax)


def get_point_epsg(lat, lon):
    """
    Get EPSG code based on latitude and longitude
    coordinates of a point
    Borrowed from geogrid.py in OPERA RTC

    Parameters
    ----------
    lat: float
        Latitude coordinate of the point
    lon: float
        Longitude coordinate of the point

    Returns
    -------
    epsg: int
        UTM zone
    """

    # "wrap" the longitude value into range [-180.0, 180.0]
    if (lon >= 180.0) or (lon <= -180.0):
        lon = (lon + 180.0) % 360.0 - 180.0

    if lat >= 75.0:
        return 3413
    elif lat <= -60.0:
        return 3031
    elif lat > 0:
        return 32601 + int(np.round((lon + 177) / 6.0))
    elif lat < 0:
        return 32701 + int(np.round((lon + 177) / 6.0))
    else:
        err_str = "'Could not determine EPSG for {0}, {1}'.format(lon, lat))"
        raise ValueError(err_str)


def make_jpl_burst_id(df):
    burst_id_jpl = (
        "t"
        + df["relative_orbit_number"].astype(str).str.zfill(3)
        + "_"
        + df["burst_id"].astype(str).str.zfill(6)
        + "_"
        + df["subswath_name"].str.lower()
    )
    return burst_id_jpl
