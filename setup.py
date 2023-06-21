"""
setup.py for OPERA burst database generator
"""

import os
import sys

from setuptools import find_packages, setup

# taken from mintpy: https://github.com/insarlab/MintPy/blob/main/setup.py
# Grab version and description from version.py
# link: https://stackoverflow.com/questions/53648900
sys.path.append(os.path.join(os.path.dirname(__file__), "src"))
from burst_db.version import release_version

LONG_DESCRIPTION = "Sentinel-1 Burst database for OPERA SAS"

setup(
    name="burst_db",
    version=release_version,
    description="Burst database for OPERA SAS",
    packages=find_packages("src"),  # include all packages under src
    package_dir={"": "src"},  # tell distutils packages are under src
    classifiers=["Programming Language :: Python"],
    url="https://github.com/opera-adt/burst_db",
    author="Seongsu Jeong; Scott J. Staniewicz",
    author_email="seongsu.jeong@jpl.nasa.gov; scott.j.staniewicz@jpl.nasa.gov",
    license=(
        "Copyright by the California Institute of Technology. ALL RIGHTS RESERVED."
    ),
    long_description=LONG_DESCRIPTION,
    # Add console scripts here
    entry_points={
        "console_scripts": [
            "opera-create-db = burst_db.make_land_frame_db:main",
        ],
    },
)
