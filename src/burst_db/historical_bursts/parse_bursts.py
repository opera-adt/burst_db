#!/usr/bin/env python
import argparse
import datetime
import logging
import os
import warnings
import zipfile
from dataclasses import astuple, dataclass
from functools import cache
from pathlib import Path
from typing import ClassVar

import boto3
import lxml.etree as ET
import numpy as np
import pandas as pd
import shapely.wkt
from dateutil.parser import parse
from eof import download
from rich.logging import RichHandler
from scipy.interpolate import InterpolatedUnivariateSpline
from shapely.geometry import LinearRing, MultiPolygon, Polygon
from tqdm.contrib.concurrent import thread_map

logger = logging.getLogger("burst_db")
h = RichHandler(rich_tracebacks=True, log_time_format="[%Y-%m-%d %H:%M:%S]")
logger.addHandler(h)
logger.setLevel(logging.INFO)

logging.getLogger("asfsmd").setLevel(logging.WARNING)


ASCENDING_NODE_TIME_TOLERANCE_IN_SEC = 1.0
T_ORBIT = (12 * 86400.0) / 175.0
PADDING_SHORT = 60

T_ORBIT = (12 * 86400.0) / 175.0


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
        track_number = int(track_number[1:])
        esa_burst_id = int(esa_burst_id)
        return cls(track_number, esa_burst_id, subswath.lower())

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
    """Extract the start/end relative orbits from manifest.safe file"""
    search_term, nsmap = _get_manifest_pattern(
        tree, ["orbitReference", "relativeOrbitNumber"]
    )
    elem_start, elem_end = tree.findall(search_term, nsmap)
    return int(elem_start.text), int(elem_end.text)


def _bursts_from_xml(annotation_path: str, orbit_path: str, open_method=open):
    """Parse bursts in Sentinel-1 annotation XML.

    Parameters:
    -----------
    annotation_path : str
        Path to Sentinel-1 annotation XML file of specific subswath and
        polarization.
    orbit_path : str
        Path the orbit file.
    open_method : function
        Function used to open annotation file.

    Returns:
    --------
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

        # Calculate ascending node crossing time from orbit;
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
def get_osv_list_from_orbit(orbit_file: str | list):
    """
    Get the list of orbit state vectors as ElementTree (ET) objects from the orbit file `orbit_file`

    - When `orbit_file` is a list, then all of the OSV lists are
      merged and sorted. The merging is done by concatenating
      the OSVs to the "base OSV list".
    - The "base OSV list" is the one that covers
      the swath's start / stop time with padding applied,
      so that the equal time spacing is guaranteed during that time period.

    Parameters
    ----------
    orbit_file: str | list
        Orbit file name, or list of the orbit file names
    swath_start: datetime.datetime
        Sensing start time of the swath
    swath_stop: datetime.datetime
        Sensing stop time of the swath

    Returns
    -------
    orbit_state_vector_list: ET
        Orbit state vector list
    """
    if isinstance(orbit_file, str):
        orbit_tree = ET.parse(orbit_file)
        orbit_state_vector_list = orbit_tree.find("Data_Block/List_of_OSVs")
        return orbit_state_vector_list

    elif isinstance(orbit_file, list) and len(orbit_file) == 1:
        orbit_tree = ET.parse(orbit_file[0])
        orbit_state_vector_list = orbit_tree.find("Data_Block/List_of_OSVs")
        return orbit_state_vector_list

    elif isinstance(orbit_file, list) and len(orbit_file) > 1:
        # Concatenate the orbit files' OSV lists
        # Assuming that the orbit was sorted, and the last orbit in the list if the one that
        # covers the sensing start / stop with padding (i.e. the base orbit file in this function)

        # Load the base orbit file
        base_orbit_tree = ET.parse(orbit_file[-1])
        orbit_state_vector_list = base_orbit_tree.find("Data_Block/List_of_OSVs")

        for i_orbit, src_orbit_filename in enumerate(orbit_file):
            if i_orbit == len(orbit_file) - 1:
                continue
            src_orbit_tree = ET.parse(src_orbit_filename)
            src_osv_vector_list = src_orbit_tree.find("Data_Block/List_of_OSVs")
            orbit_state_vector_list = merge_osv_list(
                orbit_state_vector_list, src_osv_vector_list
            )

        return orbit_state_vector_list

    else:
        raise RuntimeError(
            "Invalid orbit information provided. "
            "It has to be filename or the list of filenames."
        )


def _sort_list_of_osv(list_of_osvs):
    """
    Sort the orbit state vectors (OSV) with respect to the UTC time

    Parameters
    ----------
    list_of_osvs: ET.ElementTree
        Orbit state vectors (OSVs) as XML ElementTree (ET) objects

    Returns
    -------
    list_of_osvs: ET.ElementTree
        Sorted orbit state vectors (OSVs) with respect to UTC time
    """
    utc_osv_list = [_get_utc_time_from_osv(osv) for osv in list_of_osvs]

    sorted_index_list = [
        index for index, _ in sorted(enumerate(utc_osv_list), key=lambda x: x[1])
    ]

    list_of_osvs_copy = list_of_osvs.__copy__()

    for i_osv, _ in enumerate(list_of_osvs_copy):
        index_to_replace = sorted_index_list[i_osv]
        list_of_osvs[i_osv] = list_of_osvs_copy[index_to_replace].__copy__()

    return list_of_osvs


def merge_osv_list(list_of_osvs1, list_of_osvs2):
    """
    Merge the two orbit state vector list (OSV list) in orbit files into one.
    Apply sorting to make sure the OSVs are in chronological order
    `list_of_osvs1` will be the "base OSV list" while the OSVs in
    `list_of_osvs2` not in the time range of the base OSV list will be
    appended.

    Parameters
    ---------
        list_of_osv1: ElementTree
            Based OSV list
        list_of_osve2: ElementTree
            OSV list to be augmented to the base OSV list

    Returns
    -------
        list_of_osvs1: ElementTree
            Merged OSV list
    """
    # Extract the UTC from the OSV of the first XML
    osv1_utc_list = [_get_utc_time_from_osv(osv1) for osv1 in list_of_osvs1]

    min_utc_osv1 = min(osv1_utc_list)
    max_utc_osv1 = max(osv1_utc_list)

    for osv in list_of_osvs2.findall("OSV"):
        utc_osv2 = _get_utc_time_from_osv(osv)

        if min_utc_osv1 < utc_osv2 < max_utc_osv1:
            continue
        list_of_osvs1.append(osv)

    # sort the OSVs in the conatenated OSV
    list_of_osvs1 = _sort_list_of_osv(list_of_osvs1)

    # Adjust the count attribute in <List_of_OSVs>
    new_count = len(list_of_osvs1.findall("OSV"))
    list_of_osvs1.set("count", str(new_count))

    return list_of_osvs1


def _get_utc_time_from_osv(osv):
    """
    Extract the UTC time from orbit state vector element in orbit file

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

    Parameters:
    -----------
    tree : Element
        Element containing geolocation grid points.

    Returns:
    --------
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
    anx_time_annotation: datetime.datetime = None,
    search_length=None,
):
    """
    Estimate the time of ascending node crossing from orbit

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

    # Load the OSVs
    # utc_vec_all = []
    # pos_z_vec_all = []
    # for osv in orbit_state_vector_list:
    #     utc_str = osv.find('UTC').text.split('=')[1]
    #     utc_vec_all.append(datetime.datetime.fromisoformat(utc_str))
    #     pos_z_vec_all.append(float(osv.find('Z').text))

    # # utc_vec_all = [datetime.datetime.fromisoformat(osv.find('UTC').text.replace('UTC=',''))
    # #                for osv in orbit_state_vector_list]
    # utc_vec_all = np.array(utc_vec_all)
    # # pos_z_vec_all = [float(osv.find('Z').text)
    # #                  for osv in orbit_state_vector_list]
    # pos_z_vec_all = np.array(pos_z_vec_all)
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
        # Check if ascending node crossing has taken place;
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
            if not any([k > 180 for k in x]):
                continue

            # Otherwise, wrap longitude values down to 360 deg
            x_wrapped_minus_360 = np.asarray(x) - 360
            polys[polygon_count] = Polygon(zip(x_wrapped_minus_360, y))

        assert len(polys) == 2
    else:
        # If dateline is not crossed, treat input poly as list
        polys = [poly]

    return polys


def bursts_from_safe_dir(safe_path: str, orbit_path: str) -> list[S1Burst]:
    """Find S1Bursts in a Sentinel-1 SAFE structured directory/zipfile.

    Parameters:
    -----------
    path : str
        Path to SAFE directory.
    orbit_path : str
        Path the orbit file.

    Returns:
    --------
    bursts : list
        List of Sentinel1BurstSlc objects found in annotation XML.
    """

    def _is_zip_annotation(path: str):
        return path.split("/")[-2] == "annotation" and path.endswith(".xml")

    # find annotation file - subswath of interest
    path = Path(safe_path)
    bursts = []
    if path.suffix == ".zip":
        with zipfile.ZipFile(path, "r") as z_file:
            z_file_list = z_file.namelist()

            # find annotation file - subswath of interest
            annotation_files = [f for f in z_file_list if _is_zip_annotation(f)]

            for f_annotation in annotation_files:
                bursts.extend(
                    _bursts_from_xml(f_annotation, orbit_path, open_method=z_file.open)
                )
    elif path.is_dir():
        annotation_files = (path / "annotation").glob("s1*-iw*-slc-*")

        for f_annotation in annotation_files:
            bursts.extend(_bursts_from_xml(str(f_annotation), orbit_path))
    else:
        raise ValueError(f"Unknown file type, not a .zip or dir: {safe_path}")
    return bursts


def get_burst_rows(
    safe_file: str | Path,
    orbit_file: str | Path,
    out_dir: str | Path,
):
    try:
        outfile = (Path(out_dir) / Path(safe_file).stem).with_suffix(".csv")
        if outfile.exists():
            return
        bursts = bursts_from_safe_dir(safe_file, orbit_file)
        # all_rows.append(list(map(_to_row, bursts, safe_file))
        all_rows = [_to_row(burst, safe_file) for burst in bursts]

        pd.DataFrame(all_rows).to_csv(outfile, header=False, mode="w", index=False)
    except Exception as e:
        print(f"Failure on {safe_file}: {e}")
        outfile = (Path(out_dir) / f"failure_{Path(safe_file).stem}").touch()


def _to_row(burst: S1Burst, safe_file: str | Path):
    burst_id = burst.burst_id
    dt = burst.sensing_start.isoformat()
    # Get an approximate border, fewer points
    border = burst.border.simplify(1 / 3600).wkt
    return burst_id, dt, border, Path(safe_file).stem


PREFIX_TEMPLATE = "S1{sat}_IW_SLC__1S{sd}V_{datestr}"


def pull_safes_for_date(
    date: datetime.date,
    bucket: str,
    in_folder: str,
    satellite: str = "A",
    single_or_double: str = "D",
    out_dir: str | Path = ".",
):
    """Download SAFE files from S3 for a given date."""
    import boto3

    s3 = boto3.client("s3")
    search_term = PREFIX_TEMPLATE.format(
        sat=satellite, sd=single_or_double, datestr=date.strftime("%Y%m%d")
    )
    prefix = f"{in_folder}/{search_term}"
    print(prefix)
    response = s3.list_objects_v2(Bucket=bucket, Prefix=prefix)
    if response["KeyCount"] == 0:
        logger.info(f"No files found for {date}")
        return []
    keys = [x["Key"] for x in response["Contents"]]
    s3.download_file(bucket, keys[0], "tmp.zip")
    with zipfile.ZipFile("tmp.zip", "r") as z:
        z.extractall(out_dir)
    return [Path(x).stem for x in z.namelist()]


def make_all_safe_metadata(
    *,
    safe_list: list[str | Path] | None = None,
    out_dir: str | Path,
    orbit_file: str | Path,
    max_workers=20,
):
    print(f"Writing CSVs to {out_dir}")

    warnings.filterwarnings("ignore", category=UserWarning)  # s1reader is chatty

    def _run(file):
        get_burst_rows(file, orbit_file=orbit_file, out_dir=out_dir)

    # Use thread_map
    thread_map(
        _run,
        safe_list,
        max_workers=1,
        desc="Extracting Burst Metadata",
    )


def main() -> None:
    """Download Sentinel-1 metadata from a WKT file."""

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
        "--safe-list", help="Text file with list of SAFE granule names.", type=Path
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
        default="D",
        type=str,
        choices=["S", "D"],
        help="Single or dual polarization to download.",
    )

    args = parser.parse_args()

    # Create the output directory
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    # Get range of dates to download
    start_date = args.start_date
    end_date = args.end_date
    date_range = pd.date_range(start_date, end_date)

    # orbit_file = args.orbit_file
    # orbit_file = (
    #     "S1A_OPER_AUX_POEORB_OPOD_20210305T182757_V20150202T225944_20150204T005944.EOF"
    # )
    # make_all_safe_metadata(
    #     safe_list=Path(args.safe_list).read_text().splitlines(),
    #     out_dir=out_dir,
    #     orbit_file=orbit_file,
    # )
    # return

    # For each date, download the SAFE files
    for date in date_range:
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
        )

        logger.info("Finding bursts in SAFE files")
        make_all_safe_metadata(
            safe_list=safes,
            out_dir=out_dir,
            orbit_file=orbit_file,
            max_workers=args.max_workers,
        )   

        # Combine all the CSVs into one per date
        logger.info("Combining CSVs")
        all_csvs = list(out_dir.glob("*.csv"))
        all_rows = []
        for csv in all_csvs:
            all_rows.extend(pd.read_csv(csv, header=None).values.tolist())
            csv.unlink()
        date_output = out_dir / f"{date.strftime('%Y%m%d')}.csv"
        pd.DataFrame(all_rows).to_csv(
            date_output,
            header=False,
            mode="w",
            index=False,
        )

        # Upload to S3
        logger.info("Uploading to S3")
        s3 = boto3.client("s3")
        s3.upload_file(
            str(date_output),
            args.bucket,
            f"{args.out_folder}/{date_output.name}",
        )
        date_output.unlink()



if __name__ == "__main__":
    main()
