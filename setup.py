'''
setup.py for OPERA burst database generator
'''

import os

from setuptools import setup

__version__ = VERSION = '0.1.0'

LONG_DESCRIPTION = 'Burst database for OPERA SAS'

package_data_dict = {}

package_data_dict['rtc'] = [
    os.path.join('defaults', 'rtc_s1.yaml'),
    os.path.join('schemas', 'rtc_s1.yaml')]

setup(
    name = 'burst_db',
    version = VERSION,
    description = 'Burst database for OPERA SAS',
    package_dir = {'burst_db': 'src/burst_db'},
    include_package_data = True,
    package_data = package_data_dict,
    classifiers = ['Programming Language :: Python'],
    #scripts = ['app/rtc_s1.py'],
    install_requires = ['argparse', 'numpy', 'gdal'],
    url = 'https://github.com/opera-adt/burst_db',
    author = ('Seongsu Jeong'),
    author_email = ('seongsu.jeong@jpl.nasa.gov'),
    license = ('Copyright by the California Institute of Technology.'
               ' ALL RIGHTS RESERVED.'),
    long_description=LONG_DESCRIPTION
)
