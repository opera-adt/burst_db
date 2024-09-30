#!/usr/bin/env python
from __future__ import annotations

import argparse
import datetime
import json
import logging
import os
import shutil
import subprocess
import warnings
import zipfile
from dataclasses import astuple, dataclass
from functools import cache
from pathlib import Path
from typing import ClassVar, Optional, Sequence

import boto3
import lxml.etree as ET
import numpy as np

# import pandas as pd
import shapely.ops
import shapely.wkt
from botocore.config import Config
from botocore.exceptions import ClientError
from dateutil.parser import parse
from eof import download
from scipy.interpolate import InterpolatedUnivariateSpline
from shapely.geometry import LinearRing, MultiPolygon, Polygon
from tqdm.contrib.concurrent import thread_map

logger = logging.getLogger("burst_db")
# Make a logger good for AWS logs
h = logging.StreamHandler()
h.setFormatter(
    logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S%z",
    )
)
logger.addHandler(h)
logger.setLevel(logging.INFO)

logging.getLogger("asfsmd").setLevel(logging.WARNING)


ASCENDING_NODE_TIME_TOLERANCE_IN_SEC = 1.0
T_ORBIT = (12 * 86400.0) / 175.0
PADDING_SHORT = 60

T_ORBIT = (12 * 86400.0) / 175.0

# Use semicolon for output CSVs due to WKT commas
DELIMITER = ";"


@dataclass(frozen=True)
class S1BurstId:
    # Constants in Table 9-7 of Sentinel-1 SLC Detailed Algorithm Definition
    T_beam: ClassVar[float] = 2.758273  # interval of one burst [s]
    T_pre: ClassVar[float] = 2.299849  # Preamble time interval [s]
    T_orb: ClassVar[float] = T_ORBIT  # Nominal orbit period [s]
    track_number: int
    esa_burst_id: int
    subswath: str

    @classmethod
    def from_burst_params(
        cls,
        sensing_time: datetime.datetime,
        ascending_node_dt: datetime.datetime,
        start_track: int,
        end_track: int,
        subswath: str,
    ):
        """Calculate the unique burst ID (track, ESA burstId, swath) of a burst.

        Accounts for equator crossing frames, where the current track number of
        a burst may change mid-frame. Uses the ESA convention defined in the
        Sentinel-1 Level 1 Detailed Algorithm Definition.

        Parameters
        ----------
        sensing_time : datetime
            Sensing time of the first input line of this burst [UTC]
            The XML tag is sensingTime in the annotation file.
        ascending_node_dt : datetime
            Time of the ascending node prior to the start of the scene.
        start_track : int
            Relative orbit number at the start of the acquisition, from 1-175.
        end_track : int
            Relative orbit number at the end of the acquisition.
        subswath : str, {'IW1', 'IW2', 'IW3'}
            Name of the subswath of the burst (not case sensitive).

        Returns
        -------
        S1BurstId
            The burst ID object containing track number + ESA's burstId number + swath ID.

        Notes
        -----
        The `start_track` and `end_track` parameters are used to determine if the
        scene crosses the equator. They are the same if the frame does not cross
        the equator.

        References
        ----------
        ESA Sentinel-1 Level 1 Detailed Algorithm Definition
        https://sentinels.copernicus.eu/documents/247904/1877131/S1-TN-MDA-52-7445_Sentinel-1+Level+1+Detailed+Algorithm+Definition_v2-4.pdf/83624863-6429-cfb8-2371-5c5ca82907b8

        """
        swath_num = int(subswath[-1])
        # Since we only have access to the current subswath, we need to use the
        # burst-to-burst times to figure out
        #   1. if IW1 crossed the equator, and
        #   2. The mid-burst sensing time for IW2
        # IW1 -> IW2 takes ~0.83220 seconds
        # IW2 -> IW3 takes ~1.07803 seconds
        # IW3 -> IW1 takes ~0.84803 seconds
        burst_times = np.array([0.832, 1.078, 0.848])
        iw1_start_offsets = [
            0,
            -burst_times[0],
            -burst_times[0] - burst_times[1],
        ]
        offset = iw1_start_offsets[swath_num - 1]
        start_iw1 = sensing_time + datetime.timedelta(seconds=offset)

        start_iw1_to_mid_iw2 = burst_times[0] + burst_times[1] / 2
        mid_iw2 = start_iw1 + datetime.timedelta(seconds=start_iw1_to_mid_iw2)

        has_anx_crossing = (end_track == start_track + 1) or (
            end_track == 1 and start_track == 175
        )

        time_since_anx_iw1 = (start_iw1 - ascending_node_dt).total_seconds()
        time_since_anx = (mid_iw2 - ascending_node_dt).total_seconds()

        if (time_since_anx_iw1 - cls.T_orb) < 0:
            # Less than a full orbit has passed
            track_number = start_track
        else:
            track_number = end_track
            # Additional check for scenes which have a given ascending node
            # that's more than 1 orbit in the past
            if not has_anx_crossing:
                time_since_anx = time_since_anx - cls.T_orb

        # Eq. 9-89: ∆tb = tb − t_anx + (r - 1)T_orb
        # tb: mid-burst sensing time (sensing_time)
        # t_anx: ascending node time (ascending_node_dt)
        # r: relative orbit number   (relative_orbit_start)
        dt_b = time_since_anx + (start_track - 1) * cls.T_orb

        # Eq. 9-91 :   1 + floor((∆tb − T_pre) / T_beam )
        esa_burst_id = 1 + int(np.floor((dt_b - cls.T_pre) / cls.T_beam))

        return cls(track_number, esa_burst_id, subswath)

    @classmethod
    def from_str(cls, burst_id_str: str):
        """Parse a S1BurstId object from a string.

        Parameters
        ----------
        burst_id_str : str
            The burst ID string, e.g. "t123_000456_iw1"

        Returns
        -------
        S1BurstId
            The burst ID object containing track number + ESA's burstId number + swath ID.

        """
        track_number, esa_burst_id, subswath = burst_id_str.split("_")
        return cls(int(track_number[1:]), int(esa_burst_id), subswath.lower())

    def __str__(self):
        # Form the unique JPL ID by combining track/burst/swath
        return (
            f"t{self.track_number:03d}_{self.esa_burst_id:06d}_{self.subswath.lower()}"
        )

    def __eq__(self, other) -> bool:
        # Allows for comparison with strings, as well as S1BurstId objects
        # e.g., you can filter down burst IDs with:
        # burst_ids = ["t012_024518_iw3", "t012_024519_iw3"]
        # bursts = [b for b in bursts if b.burst_id in burst_ids]
        if isinstance(other, str):
            return str(self) == other
        elif isinstance(other, S1BurstId):
            return astuple(self) == astuple(other)
        else:
            raise ValueError(f"Cannot compare to {other}")


@dataclass(frozen=True)
class S1Burst:
    burst_id: S1BurstId
    sensing_start: datetime.datetime
    border: MultiPolygon


def as_datetime(t_str: str) -> datetime.datetime:
    return datetime.datetime.fromisoformat(t_str)


def _get_manifest_pattern(tree: ET, keys: list):
    """Get the search path to extract data from the ET of manifest.safe.

    Parameters
    ----------
    tree : xml.etree.ElementTree
        ElementTree containing the parsed 'manifest.safe' file.
    keys : list
        List of keys to search for in the manifest file.

    Returns
    -------
    str
        Search path to extract from the ET of the manifest.safe XML.

    """
    # https://lxml.de/tutorial.html#namespaces
    # Get the namespace from the root element to avoid full urls
    # This is a dictionary with a short name containing the labels in the
    # XML tree, and the fully qualified URL as the value.
    try:
        nsmap = tree.nsmap
    except AttributeError:
        nsmap = tree.getroot().nsmap
    # path to xmlData in manifest
    xml_meta_path = "metadataSection/metadataObject/metadataWrap/xmlData"
    safe_terms = "/".join([f"safe:{key}" for key in keys])
    return f"{xml_meta_path}/{safe_terms}", nsmap


def get_start_end_track(tree: ET):
    """Extract the start/end relative orbits from manifest.safe file."""
    search_term, nsmap = _get_manifest_pattern(
        tree, ["orbitReference", "relativeOrbitNumber"]
    )
    elem_start, elem_end = tree.findall(search_term, nsmap)
    return int(elem_start.text), int(elem_end.text)


def _bursts_from_xml(annotation_path: str, orbit_path: str, open_method=open):
    """Parse bursts in Sentinel-1 annotation XML.

    Parameters
    ----------
    annotation_path : str
        Path to Sentinel-1 annotation XML file of specific subswath and
        polarization.
    orbit_path : str
        Path the orbit file.
    open_method : function
        Function used to open annotation file.

    Returns
    -------
    bursts : list
        List of Sentinel1BurstSlc objects found in annotation XML.

    """
    _, tail = os.path.split(annotation_path)
    platform_id, subswath, _, pol = [x.upper() for x in tail.split("-")[:4]]

    # parse manifest.safe to retrieve IPF version
    # Also load the Product annotation - for EAP calibration and RFI information
    manifest_path = (
        os.path.dirname(annotation_path).replace("annotation", "") + "manifest.safe"
    )
    tree_manifest = ET.parse(manifest_path)
    # Parse out the start/end track to determine if we have an
    # equator crossing (for the burst_id calculation).
    start_track, end_track = get_start_end_track(tree_manifest)

    # Nearly all metadata loaded here is common to all bursts in annotation XML
    tree = ET.parse(annotation_path)

    image_info_element = tree.find("imageAnnotation/imageInformation")
    ascending_node_time_annotation = as_datetime(
        image_info_element.find("ascendingNodeTime").text
    )
    first_line_utc_time = as_datetime(
        image_info_element.find("productFirstLineUtcTime").text
    )

    # orbit_number = int(tree.find('adsHeader/absoluteOrbitNumber').text)
    boundary_polys = _get_burst_bounds(tree)

    # find orbit state vectors in 'Data_Block/List_of_OSVs'
    if orbit_path:
        orbit_state_vector_list = get_osv_list_from_orbit(orbit_path)

        # Calculate ascending node crossing time from orbit,
        # compare with the info from annotation
        try:
            ascending_node_time_orbit = get_ascending_node_time_orbit(
                orbit_state_vector_list,
                first_line_utc_time,
                ascending_node_time_annotation,
            )
            ascending_node_time = ascending_node_time_orbit

        except ValueError:
            warnings.warn(
                "Cannot estimate ascending node crossing time from Orbit. "
                "Using the ascending node time in annotation."
            )
            ascending_node_time_orbit = None
            ascending_node_time = ascending_node_time_annotation

        # Calculate the difference in the two ANX times, and give
        # warning message when the difference is noticeable.
        if ascending_node_time_orbit is not None:
            diff_ascending_node_time_seconds = (
                ascending_node_time_orbit - ascending_node_time_annotation
            ).total_seconds()
            if (
                abs(diff_ascending_node_time_seconds)
                > ASCENDING_NODE_TIME_TOLERANCE_IN_SEC
            ):
                warnings.warn(
                    "ascending node time error larger than "
                    f"{ASCENDING_NODE_TIME_TOLERANCE_IN_SEC} seconds was detected: "
                    f"Error = {diff_ascending_node_time_seconds} seconds."
                )

    else:
        warnings.warn(
            "Orbit file was not provided. "
            "Using the ascending node time from annotation."
        )
        ascending_node_time = ascending_node_time_annotation
        orbit_state_vector_list = []

    # load individual burst
    burst_list_elements = tree.find("swathTiming/burstList")

    bursts = []
    for burst_list_element, poly in zip(burst_list_elements, boundary_polys):
        # Zero Doppler azimuth time of the first line of this burst
        sensing_start = as_datetime(burst_list_element.find("azimuthTime").text)

        # # Sensing time of the first input line of this burst [UTC]
        sensing_time = as_datetime(burst_list_element.find("sensingTime").text)
        # Create the burst ID to match the ESA ID scheme
        burst_id = S1BurstId.from_burst_params(
            sensing_time, ascending_node_time, start_track, end_track, subswath
        )
        bursts.append(
            S1Burst(
                burst_id=burst_id,
                sensing_start=sensing_start,
                border=poly,
            )
        )

    return bursts


@cache
def get_osv_list_from_orbit(orbit_file: str):
    """Get the orbit state vectors as ET objects from `orbit_file`.

    Parameters
    ----------
    orbit_file: str | list
        Orbit file name, or list of the orbit file names

    Returns
    -------
    orbit_state_vector_list: ET
        Orbit state vector list

    """
    orbit_tree = ET.parse(orbit_file)
    orbit_state_vector_list = orbit_tree.find("Data_Block/List_of_OSVs")
    return orbit_state_vector_list


def _get_utc_time_from_osv(osv):
    """Extract the UTC time from orbit state vector element in orbit file.

    Parameters
    ----------
    osv: ElementTree
        orbit state vector parsed as .xml

    Returns
    -------
    datetime_utc: datetime.datetime
        Orbit state vector's UTC time

    """
    utc_osv_string = osv.find("UTC").text.replace("UTC=", "")
    return as_datetime(utc_osv_string)


def _get_burst_bounds(tree):
    """Parse grid points list and get border.

    Parameters
    ----------
    tree : Element
        Element containing geolocation grid points.

    Returns
    -------
    boundary_polys : list[MultiPolygon]
        List of burst boundaries as shapely MultiPolygons

    """
    # find element tree
    grid_pt_list = tree.find("geolocationGrid/geolocationGridPointList")

    # read in all points
    n_grid_pts = int(grid_pt_list.attrib["count"])
    lines = np.empty(n_grid_pts)
    pixels = np.empty(n_grid_pts)
    lats = np.empty(n_grid_pts)
    lons = np.empty(n_grid_pts)
    for i, grid_pt in enumerate(grid_pt_list):
        lines[i] = int(grid_pt[2].text)
        pixels[i] = int(grid_pt[3].text)
        lats[i] = float(grid_pt[4].text)
        lons[i] = float(grid_pt[5].text)

    unique_line_indices = np.unique(lines)
    boundary_polys = []

    # zip lines numbers of bursts together and iterate
    for ln0, ln1 in zip(unique_line_indices[:-1], unique_line_indices[1:]):
        # create masks for lines in current burst
        mask0 = lines == ln0
        mask1 = lines == ln1

        # reverse order of 2nd set of points so plots of boundaries
        # are not connected by a diagonal line
        burst_lons = np.concatenate((lons[mask0], lons[mask1][::-1]))
        burst_lats = np.concatenate((lats[mask0], lats[mask1][::-1]))

        poly = Polygon(zip(burst_lons, burst_lats))
        boundary_polys.append(MultiPolygon(check_dateline(poly)))

    return boundary_polys


@cache
def _get_utc_z(orbit_state_vector_list):
    utc_vec_all = []
    pos_z_vec_all = []
    for osv in orbit_state_vector_list:
        utc_str = osv.find("UTC").text.split("=")[1]
        utc_vec_all.append(datetime.datetime.fromisoformat(utc_str))
        pos_z_vec_all.append(float(osv.find("Z").text))
    utc_vec_all = np.array(utc_vec_all)
    pos_z_vec_all = np.array(pos_z_vec_all)
    return utc_vec_all, pos_z_vec_all


def get_ascending_node_time_orbit(
    orbit_state_vector_list: ET,
    sensing_time: datetime.datetime,
    anx_time_annotation: Optional[datetime.datetime] = None,
    search_length=None,
):
    """Estimate the time of ascending node crossing from orbit.

    Parameters
    ----------
    orbit_state_vector_list: ET
        XML elements that points to the list of orbit information.
        Each element should contain the information below:
        TAI, UTC, UT1, Absolute_Orbit, X, Y, Z, VX, VY, VZ, and Quality

    sensing_time: datetime.datetime
        Sensing time of the data

    anx_time_annotation: datetime.datetime
        ascending node crossing (ANX) time retrieved from annotation data.
        - If it is provided, then the orbit ANX time close to this will be returned.
        - If it is `None`, then the orbit ANX time closest to (but no later than)
          the sensing time will be returned.

    search_length: datetime.timedelta
        Search length in time to crop the data in `orbit_state_vector_list`

    Returns
    -------
    _ : datetime.datetime
        Ascending node crossing time calculated from orbit

    """
    # Crop the orbit information before 2 * (orbit period) of sensing time
    if search_length is None:
        search_length = datetime.timedelta(seconds=2 * T_ORBIT)

    utc_vec_all, pos_z_vec_all = _get_utc_z(orbit_state_vector_list)

    # NOTE: tried to apply the same amount of pading in `get_burst_orbit` to
    # get as similar results as possible.
    pad = datetime.timedelta(seconds=PADDING_SHORT)

    # Constrain the search area
    flag_search_area = (utc_vec_all > sensing_time - search_length - pad) & (
        utc_vec_all < sensing_time + pad
    )
    orbit_time_vec = utc_vec_all[flag_search_area]
    orbit_z_vec = pos_z_vec_all[flag_search_area]

    # Detect the event of ascending node crossing from the cropped orbit info.
    # The algorithm was inspired by Scott Staniewicz's PR in the link below:
    # https://github.com/opera-adt/s1-reader/pull/120/
    datetime_ascending_node_crossing_list = []

    # Iterate through the z coordinate in orbit object to
    # detect the ascending node crossing
    iterator_z_time = zip(orbit_z_vec, orbit_z_vec[1:], orbit_time_vec)
    for index_z, (z_prev, z, _) in enumerate(iterator_z_time):
        # Check if ascending node crossing has taken place:
        if not z_prev < 0 <= z:
            continue

        # Extract z coords. and time for interpolation
        index_from = max(index_z - 3, 0)
        index_to = min(index_z + 3, len(orbit_z_vec))
        z_around_crossing = orbit_z_vec[index_from:index_to]
        time_around_crossing = orbit_time_vec[index_from:index_to]

        # Set up spline interpolator and interpolate the time when z is equal to 0.0
        datetime_orbit_ref = time_around_crossing[0]
        relative_time_around_crossing = [
            (t - datetime_orbit_ref).total_seconds() for t in time_around_crossing
        ]

        interpolator_time = InterpolatedUnivariateSpline(
            z_around_crossing, relative_time_around_crossing, k=1
        )
        relative_t_interp = interpolator_time(0.0)

        # Convert the interpolated time into datetime
        # Add the ascending node crossing time if it's before the sensing time
        datetime_ascending_crossing = datetime_orbit_ref + datetime.timedelta(
            seconds=float(relative_t_interp)
        )
        if datetime_ascending_crossing < sensing_time:
            datetime_ascending_node_crossing_list.append(datetime_ascending_crossing)

    if len(datetime_ascending_node_crossing_list) == 0:
        raise ValueError(
            "Cannot detect ascending node crossings "
            "from the orbit information provided."
        )

    # Return the most recent time for ascending node crossing
    if anx_time_annotation is None:
        anx_time_orbit = max(datetime_ascending_node_crossing_list)
    elif isinstance(anx_time_annotation, datetime.datetime):
        anx_time_orbit = min(
            datetime_ascending_node_crossing_list,
            key=lambda anxtime: abs(anxtime - anx_time_annotation),
        )
    else:
        raise ValueError(
            "datatype of `anx_time_annotation` has to be `datetime.datetime`"
        )

    return anx_time_orbit


def check_dateline(poly):
    """Split `poly` if it crosses the dateline.

    Parameters
    ----------
    poly : shapely.geometry.Polygon
        Input polygon.

    Returns
    -------
    polys : list of shapely.geometry.Polygon
         A list containing: the input polygon if it didn't cross
        the dateline, or two polygons otherwise (one on either
        side of the dateline).

    """
    xmin, _, xmax, _ = poly.bounds
    # Check dateline crossing
    if (xmax - xmin) > 180.0:
        dateline = shapely.wkt.loads("LINESTRING( 180.0 -90.0, 180.0 90.0)")

        # build new polygon with all longitudes between 0 and 360
        x, y = poly.exterior.coords.xy
        new_x = (k + (k <= 0.0) * 360 for k in x)
        new_ring = LinearRing(zip(new_x, y))

        # Split input polygon
        # (https://gis.stackexchange.com/questions/232771/splitting-polygon-by-linestring-in-geodjango_)
        merged_lines = shapely.ops.linemerge([dateline, new_ring])
        border_lines = shapely.ops.unary_union(merged_lines)
        decomp = shapely.ops.polygonize(border_lines)

        polys = list(decomp)

        # The Copernicus DEM used for NISAR processing has a longitude
        # range [-180, +180]. The current version of gdal.Translate
        # does not allow to perform dateline wrapping. Therefore, coordinates
        # above 180 need to be wrapped down to -180 to match the Copernicus
        # DEM longitude range
        for polygon_count in range(2):
            x, y = polys[polygon_count].exterior.coords.xy
            if not any(k > 180 for k in x):
                continue

            # Otherwise, wrap longitude values down to 360 deg
            x_wrapped_minus_360 = np.asarray(x) - 360
            polys[polygon_count] = Polygon(zip(x_wrapped_minus_360, y))

        assert len(polys) == 2
    else:
        # If dateline is not crossed, treat input poly as list
        polys = [poly]

    return polys


def unzip_safe(file_path: str | Path, output_directory: Path):
    # Extracts the zip file ensuring a consistent folder structure

    safe_name = Path(file_path).stem
    with zipfile.ZipFile(file_path, "r") as zip_ref:
        # If the zip starts with "safe_folders/", then adjust extraction
        # get the number of parts in the `manifest`
        manifest_parts = next(
            Path(f).parts for f in zip_ref.namelist() if "manifest.safe" in f
        )
        manifest_level = len(manifest_parts)

        # Determine extraction logic based on contents
        if manifest_level == 2:
            # Extract everything to the output directory
            zip_ref.extractall(output_directory)
        elif manifest_level == 3:
            # Extract , then after move up one directory
            zip_ref.extractall(output_directory)
            # Get the subdirectory name (file name without .zip)
            shutil.move(
                output_directory / manifest_parts[0] / manifest_parts[1],
                output_directory,
            )
            # extraction_path = output_directory / safe_name
            # # Move the contents of the subdirectory up one level
            # for f in extraction_path.iterdir():
            #     shutil.move(f, output_directory)
        else:
            # Determine the subdirectory name (file name without .zip)
            extraction_path = output_directory / safe_name
            zip_ref.extractall(extraction_path)


def bursts_from_safe_dir(
    safe_path: str | Path, orbit_path: str | Path
) -> list[S1Burst]:
    """Find S1Bursts in a Sentinel-1 SAFE structured directory/zipfile.

    Parameters
    ----------
    path : str | Path
        Path to SAFE directory.
    orbit_path : str | Path
        Path the orbit file.

    Returns
    -------
    bursts : list
        List of Sentinel1BurstSlc objects found in annotation XML.

    """

    def _is_zip_annotation(path: str):
        return "annotation" in path and path.endswith(".xml")

    # find annotation file - subswath of interest
    path = Path(safe_path)
    bursts = []
    if path.suffix == ".zip":
        unzip_safe(safe_path, path.parent)
        path = path.parent / path.stem

    annotation_files = (path / "annotation").glob("s1*-iw*-slc-*")
    if not annotation_files:
        raise ValueError(f"No annotation files found in {path}")

    for f_annotation in annotation_files:
        bursts.extend(_bursts_from_xml(str(f_annotation), str(orbit_path)))
    if not bursts:
        logger.error(f"No bursts found in {path}")

    # else:
    #     raise ValueError(f"Unknown file type, not a .zip or dir: {safe_path}")
    return bursts


def get_burst_rows(
    safe_file: str | Path,
    orbit_file: str | Path,
    out_dir: str | Path,
) -> Path | None:
    try:
        outfile = (Path(out_dir) / Path(safe_file).stem).with_suffix(".csv")
        if outfile.exists():
            return outfile
        bursts = bursts_from_safe_dir(safe_file, orbit_file)
        # all_rows.append(list(map(_to_row, bursts, safe_file))
        logger.debug(f"Found {len(bursts)} bursts in {safe_file}")
        all_rows = [_to_row(burst, safe_file) for burst in bursts]

        # pd.DataFrame(all_rows).to_csv(outfile, header=False, mode="w", index=False)
        _to_csv(all_rows, outfile)
    except Exception:
        logger.error(f"Failure on {safe_file}", exc_info=True)
        return None
        # outfile = (Path(out_dir) / f"failure_{Path(safe_file).stem}")
        # outfile.touch()
    return outfile


def _to_row(burst: S1Burst, safe_file: str | Path):
    burst_id = burst.burst_id
    dt = burst.sensing_start.isoformat()
    # Get an approximate border, fewer points
    border = burst.border.simplify(1 / 3600).wkt
    return burst_id, dt, border, Path(safe_file).stem


def _to_csv(row_list: list, outfile: str | Path):
    """Write a csv without pandas."""
    import csv

    logger.debug(f"Writing {len(row_list)} rows to {outfile}")
    with open(outfile, "w") as f:
        writer = csv.writer(f, delimiter=DELIMITER)
        writer.writerows(row_list)


def pull_safes_for_date(
    date: datetime.date,
    bucket: str,
    in_folder: str,
    satellite: str = "A",
    single_or_double: str = "*",
    out_dir: str | Path = ".",
    max_workers: int = 3,
    full_safe_list: str | Path | None = None,
) -> list[Path]:
    """Download SAFE files from S3 for a given date.

    If `full_safe_list` is passed, will search the text file
    for the matching safe names to download.
    Otherwise, uses `s5cmd` to perform a prefix search and `cp`.
    """
    search_template = "S1{sat}_IW_SLC__1S{sd}V_{datestr}"
    search_term = search_template.format(
        sat=satellite, sd=single_or_double, datestr=date.strftime("%Y%m%d")
    )
    if not full_safe_list:
        return _s5cmd_copy(search_term, bucket, in_folder, out_dir, max_workers)
    else:
        matching_granules = _find_matching(search_term, full_safe_list)
        zip_filenames = [f"{f}.SAFE.zip" for f in matching_granules]
        return _get_objects(
            zip_filenames, bucket, in_folder, Path(out_dir), max_workers=max_workers
        )


def _s5cmd_copy(
    search_term: str,
    bucket: str,
    in_folder: str,
    out_dir: str | Path = ".",
    max_workers: int = 10,
) -> list[Path]:
    prefix = f"{in_folder}/{search_term}"
    full_s3_path = f"s3://{bucket}/{prefix}"
    cmd = [
        "s5cmd",
        "--json",
        "--numworkers",
        str(max_workers),
        "cp",
        f"{full_s3_path}*",
        str(Path(out_dir).resolve()),
    ]
    out = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        check=True,
    )
    loaded = list(map(json.loads, out.stdout.splitlines()))
    if not loaded:
        logger.info(f"No files found matching {search_term}")
        return []
    loaded_paths = []
    for d in loaded:
        if not d["success"]:
            logger.error(f"Failed to download {d['source']}")
        else:
            loaded_paths.append(Path(d["destination"]))
    return loaded_paths

    # s3 = boto3.client("s3")
    # response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    # if response["KeyCount"] == 0:
    #     logger.info(f"No files found for {date}")
    #     return []
    # keys = [x["Key"] for x in response["Contents"]]
    # s3.download_file(bucket, keys[0], "tmp.zip")
    # with zipfile.ZipFile("tmp.zip", "r") as z:
    #     z.extractall(out_dir)
    # return [Path(x).stem for x in z.namelist()]


def _find_matching(search_term: str, full_safe_list: str | Path) -> list[str]:
    out = subprocess.run(
        [
            "grep",
            "-G",  # Use regex grep
            search_term.replace("*", "."),  # need '.' for regex match-all
            str(full_safe_list),
        ],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if out.returncode == 1:
        logger.error(f"Failed to grep {search_term}")
        return []
    elif out.returncode > 1:
        logger.error(f"Failed to grep {search_term}: {out.stderr}")
        out.check_returncode()
    return out.stdout.splitlines()


def _get_objects(
    zip_filenames: list[str],
    bucket: str,
    s3_folder: str,
    out_dir: Path,
    max_workers: int = 3,
) -> list[Path]:
    from itertools import repeat

    def _download(s3_filename) -> Path | None:
        s3_client, filename = s3_filename
        output_file = out_dir / filename
        key = f"{s3_folder}/{filename}"
        try:
            s3_client.download_file(Bucket=bucket, Key=key, Filename=output_file)
        except ClientError as e:
            if "Not Found" in str(e):
                logger.exception(f"s3://{bucket}/{key} doesn't exist")
            else:
                logger.error(f"Failed to download s3://{bucket}/{key}", exc_info=True)
            return None
        except Exception:
            logger.error(f"Failed to download s3://{bucket}/{key}", exc_info=True)
            return None
        return output_file

    s3 = boto3.client("s3")
    outs = thread_map(
        _download,
        zip(repeat(s3), zip_filenames),
        max_workers=max_workers,
        total=len(zip_filenames),
    )
    return [f for f in outs if f is not None]


def make_all_safe_metadata(
    *,
    safe_list: Sequence[str | Path],
    out_dir: str | Path,
    orbit_file: str | Path,
) -> list[Path]:
    warnings.filterwarnings("ignore", category=UserWarning)  # s1reader is chatty

    # def _run(file):
    #     get_burst_rows(file, orbit_file=orbit_file, out_dir=out_dir)
    # Use thread_map
    # thread_map(
    #     _run,
    #     safe_list,
    #     max_workers=1,
    #     desc="Extracting Burst Metadata",
    # )
    csv_files = []
    for file in safe_list:
        # _run(file)
        out = get_burst_rows(file, orbit_file=orbit_file, out_dir=out_dir)
        if out is not None:
            csv_files.append(out)
    return csv_files


def _combine_csvs_by_date(
    csv_files: list[Path],
    output: str | Path,
    no_clean: bool = False,
) -> None:
    # Read each of the files in `csv_files` and simply concatenate them to
    # a single file in `date_output`
    # Open the output file in write mode
    with open(output, "w") as outfile:
        # Loop over each text file and write its contents to the output file
        for p in csv_files:
            with open(p, "r") as infile:
                outfile.write(infile.read())
            if not no_clean:
                logger.debug("Removing CSVs")
                p.unlink()


def _is_valid_date(
    check_date: datetime.datetime | datetime.date, satellite: str
) -> bool:
    if isinstance(check_date, datetime.datetime):
        comparison_date = check_date.date()
    else:
        comparison_date = check_date

    # POEORB files get generated about 3 weeks after the acquisition date
    # TODO: do i care? should i allow resorbs?
    max_date = (datetime.datetime.today() - datetime.timedelta(days=21)).date()
    # max_date = datetime.date(2023, 8, 1)
    # https://cmr.earthdata.nasa.gov/cloudstac/ASF/collections/SENTINEL-1A_SLC.v1/2014/10
    s1a_start_date = datetime.date(2014, 10, 3)
    # https://cmr.earthdata.nasa.gov/cloudstac/ASF/collections/SENTINEL-1B_SLC.v1/2016/08
    s1b_start_date = datetime.date(2016, 8, 20)
    # https://cmr.earthdata.nasa.gov/cloudstac/ASF/collections/SENTINEL-1B_SLC.v1/2021/12
    s1b_end_date = datetime.date(2021, 12, 23)
    if comparison_date > max_date:
        logger.info(f"Skipping {comparison_date} since it is past {max_date = }")
        return False
    if satellite == "A":
        if comparison_date < s1a_start_date:
            logger.info(
                f"Skipping {comparison_date} since it is before {s1a_start_date = }"
            )
            return False
    else:
        if comparison_date < s1b_start_date:
            logger.info(
                f"Skipping {comparison_date} since it is before {s1b_start_date = }"
            )
            return False
        if comparison_date > s1b_end_date:
            logger.info(
                f"Skipping {comparison_date} since it is after {s1b_end_date = }"
            )
            return False
    return True


def _get_s3_key(
    file_to_upload, out_folder, satellite, date, combine_by_date=True
) -> str:
    if combine_by_date:
        # Only nest by year
        year = date.strftime("%Y")
        key = f"{out_folder}/S1{satellite}/{year}/{file_to_upload.name}"
    else:
        # title by satellite / date
        year_month = date.strftime("%Y/%m")
        key = f"{out_folder}/S1{satellite}/{year_month}/{file_to_upload.name}"
    return key


def _get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download S1 metadata from a list of SAFE granule names.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--out-dir", default="scratch", type=Path, help="Temporary output directory."
    )
    parser.add_argument("--bucket", help="S3 bucket name to upload results.")
    parser.add_argument(
        "--in-folder", help="S3 folder name within bucket where inputs are located."
    )
    parser.add_argument(
        "--out-folder", help="S3 folder name within bucket to upload results."
    )
    parser.add_argument(
        "--start-date",
        default="2014-10-01",
        type=lambda x: parse(x),
        help="Date to start downloading SAFE files.",
    )
    parser.add_argument(
        "--end-date",
        default=datetime.datetime.today(),
        type=lambda x: parse(x),
        help="Date to stop downloading SAFE files.",
    )
    parser.add_argument(
        "--max-workers",
        default=3,
        type=int,
        help=(
            "Number of workers to use for downloading SAFE metadata. "
            "Each worker will get on of the batches of SAFE granules."
        ),
    )
    parser.add_argument(
        "--satellite",
        default="A",
        type=str,
        choices=["A", "B"],
        help="Sentinel-1 satellite to download.",
    )
    parser.add_argument(
        "--single-or-double",
        default="*",
        type=str,
        choices=["*", "S", "D"],
        help="Single or dual polarization to download. '*' pulls both.",
    )
    parser.add_argument(
        "--no-clean",
        action="store_true",
        help="Don't clean up SAFE files after processing.",
    )
    parser.add_argument(
        "--combine-by-date",
        action="store_true",
        help=(
            "Concatenate separate SAFE csvs into one per date. "
            "Saves file named YYYYMMDD.csv"
        ),
    )
    parser.add_argument(
        "--full-safe-list",
        type=Path,
        help=(
            "Path to a text file containing the full list of SAFE granules. "
            "If provided, will use this to search for SAFE files instead of "
            "using `s5cmd` (and is much faster)."
        ),
    )
    parser.add_argument(
        "--skip-if-exists",
        action="store_true",
        help="Don't reprocess the YYYYMMDD.csv file if it's on S3 already.",
    )
    parser.add_argument(
        "--no-upload", action="store_true", help="Don't upload the final CSVs to S3"
    )
    return parser


def main() -> None:
    """Download Sentinel-1 metadata from a WKT file."""
    parser = _get_parser()
    args = parser.parse_args()

    # Create the output directory
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Get range of dates to download
    start_date = args.start_date
    end_date = args.end_date
    # date_range = pd.date_range(start_date, end_date)
    date = start_date
    date_range = [date]
    while date < end_date:
        date += datetime.timedelta(days=1)
        date_range.append(date)

    config = Config(retries={"max_attempts": 10, "mode": "standard"})
    s3 = boto3.client("s3", config=config)

    # TEMP: remove later. unzip the file if it exists
    if not args.full_safe_list and Path("all_cmr_stac_safes_sorted.txt.zip").exists():
        subprocess.run(["unzip", "all_cmr_stac_safes_sorted.txt.zip"], check=False)
        args.full_safe_list = Path("all_cmr_stac_safes_sorted.txt")

    # For each date, download the SAFE files
    for date in date_range:
        if not _is_valid_date(date, args.satellite):
            continue

        # Check for the current date key, skip if it exists
        if args.combine_by_date and args.skip_if_exists:
            key = _get_s3_key(
                Path(f"{date.strftime('%Y%m%d')}.csv.zip"),
                args.out_folder,
                args.satellite,
                date,
                combine_by_date=args.combine_by_date,
            )
            try:
                s3.head_object(Bucket=args.bucket, Key=key)
                logger.info(f"Skipping {date}, {key} already exists")
                continue
            except ClientError as e:
                if e.response["Error"]["Message"] != "Not Found":
                    raise
                # Otherwise, continue. Expected error is 404 not found

        logger.info(
            f"Pulling SAFE folder for {date}, S1{args.satellite},"
            f" S{args.single_or_double}V pol"
        )
        safes = pull_safes_for_date(
            date=date,
            bucket=args.bucket,
            in_folder=args.in_folder,
            out_dir=out_dir,
            satellite=args.satellite,
            single_or_double=args.single_or_double,
            full_safe_list=args.full_safe_list,
        )
        logger.info(
            f"Found {len(safes)} SAFE files for {date}, S1{args.satellite},"
            f" S{args.single_or_double}V pol"
        )
        if not safes:
            logger.info(f"No SAFE files found for {date}, skipping")
            continue

        logger.info(f"Downloading orbit file for {date}")
        orbit_file = download.main(
            save_dir=out_dir,
            date=date_range[0],
            mission=f"S1{args.satellite}",
        )[0]

        logger.info("Finding bursts in SAFE files")
        csv_files = make_all_safe_metadata(
            safe_list=safes,
            out_dir=out_dir,
            orbit_file=orbit_file,
            # max_workers=args.max_workers,
        )

        # Combine all the CSVs into one per date
        if args.combine_by_date:
            date_output = out_dir / f"{date.strftime('%Y%m%d')}.csv"
            logger.info("Combining CSVs into %s", date_output)
            _combine_csvs_by_date(csv_files, output=date_output, no_clean=args.no_clean)
            # We only need to upload the combined CSV
            csv_files = [date_output]

        if args.no_upload:
            logger.info("Skipping upload to S3.")
            return

        logger.info(f"Uploading {len(csv_files)} to S3")
        for file in csv_files:
            logger.debug("Zipping output")
            # file_to_upload = date_output.with_suffix(".csv.zip")
            file_to_upload = file.with_suffix(".csv.zip")
            with zipfile.ZipFile(file_to_upload, "w") as z:
                # arcname is the name of the file inside the zip
                z.write(file, arcname=file.name)
            file.unlink()

            key = _get_s3_key(file_to_upload, args.out_folder, args.satellite, date)
            logger.info(f"Uploading {key} to s3://{args.bucket}")
            s3.upload_file(str(file_to_upload), args.bucket, key)

        if not args.no_clean:
            for file in csv_files:
                file.with_suffix(".csv.zip").unlink()
            # Clean up the SAFE files
            Path(orbit_file).unlink()
            for file in out_dir.glob("*.SAFE*"):
                if file.is_dir():
                    shutil.rmtree(file)
                else:
                    file.unlink()


if __name__ == "__main__":
    main()
