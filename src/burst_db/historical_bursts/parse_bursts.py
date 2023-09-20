#!/usr/bin/env python
import datetime
import os
import warnings
import zipfile
from dataclasses import dataclass
from pathlib import Path

import lxml.etree as ET
import numpy as np
from nisar.workflows.stage_dem import check_dateline
from shapely.geometry import MultiPolygon, Polygon

from s1reader.s1_burst_id import S1BurstId
from s1reader.s1_reader import (
    ASCENDING_NODE_TIME_TOLERANCE_IN_SEC,
    as_datetime,
    get_ascending_node_time_orbit,
    get_osv_list_from_orbit,
    get_start_end_track,
)


@dataclass(frozen=True)
class S1Burst:
    burst_id: S1BurstId
    sensing_start: datetime.datetime
    border: MultiPolygon


def burst_id_from_xml(annotation_path: str, orbit_path: str, open_method=open):
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
    boundary_polys = get_burst_boundaries(tree)

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
            annotation_files = [
                f for f in z_file_list if _is_zip_annotation(f)
            ]

            for f_annotation in annotation_files:
                bursts.extend(
                    burst_id_from_xml(f_annotation, orbit_path, open_method=z_file.open)
                )
    elif path.is_dir():
        annotation_files = (path / "annotation").glob("s1*-iw*-slc-*")

        for f_annotation in annotation_files:
            bursts.extend(burst_id_from_xml(str(f_annotation), orbit_path))
    else:
        raise ValueError(f"Unknown file type, not a .zip or dir: {safe_path}")
    return bursts


def get_burst_boundaries(tree):
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
    for (ln0, ln1) in zip(unique_line_indices[:-1], unique_line_indices[1:]):
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