import subprocess
from pathlib import Path

from burst_db import VERSION_CLEAN, __version__


def create_2d_geojsons(
    in_file: str = f"opera-s1-disp-{VERSION_CLEAN}.gpkg",
) -> list[Path]:
    """Create 2D GeoJSON files from a 3D GeoPackage database.

    This function takes a 3D GeoPackage database containing burst geometry
    and metadata, extracts 2D representations of the geometry.
    The 3rd Z dimension is not used and all 0s from the ESA database,
    and causes some issues with the GeoJSON files.

    OGR commands are used to convert the 3D data to 2D. SQL queries filter
    for land bursts and simplify the geometries.

    GeoJSON files are output for frame and burst geometries.
    The files are zipped up after creation.
    """
    # Get the burst_db version
    print(f"Burst database software version: {__version__}")

    db_file_2d = f"opera-s1-disp-{VERSION_CLEAN}-2d.gpkg"
    # Check if db_file_2d exists
    if not Path(db_file_2d).is_file():
        print(f"Creating {db_file_2d}")
        query_base = "SELECT AsGPB(CastToXY(GeomFromGPB(geom))) AS geom, * from"
        ogr2ogr_base_cmd = "ogr2ogr -nln frames -nlt MULTIPOLYGON -dialect sqlite"

        # Run ogr2ogr commands
        subprocess.run(
            f'{ogr2ogr_base_cmd} -sql "{query_base} frames" "{db_file_2d}" "{in_file}"',
            shell=True,
            check=True,
        )
        for table_name in ["burst_id_map", "frames_bursts"]:
            subprocess.run(
                f'ogr2ogr -append -nln {table_name} -dialect sqlite -sql "SELECT * from'
                f' {table_name}" "{db_file_2d}" "{in_file}"',
                shell=True,
                check=True,
            )

    # SQL queries
    sql_query_burst_id = """
SELECT b.fid,
    b.burst_id_jpl,
    is_land,
    is_north_america,
    b.orbit_pass,
    st_simplify(b.geom, 0.1) AS geometry
FROM frames f
    JOIN frames_bursts fb ON fb.frame_fid = f.fid
    JOIN burst_id_map b ON fb.burst_ogc_fid = b.fid
;
"""

    sql_query_frame = """
SELECT f.fid,
    is_land,
    is_north_america,
    orbit_pass,
    st_simplify(f.geom, 0.1) AS geometry
FROM frames f
;
"""

    output_files = []
    # Run ogr2ogr for GeoJSON and zip the files
    for query, name in [
        (sql_query_frame, "frame-geometries-simple"),
        (sql_query_burst_id, "burst-id-geometries-simple"),
    ]:
        json_file = f"{name}-{VERSION_CLEAN}.geojson"
        zip_file = f"{json_file}.zip"
        cmd = (
            f'ogr2ogr -f GeoJSON -preserve_fid -dialect sqlite -sql "{query}"'
            f' "{json_file}" {db_file_2d}'
        )
        print(cmd)
        subprocess.run(
            cmd,
            shell=True,
            check=True,
        )
        subprocess.run(f"zip {zip_file} {json_file}", shell=True, check=True)
        output_files.append(Path(zip_file))

    return output_files
