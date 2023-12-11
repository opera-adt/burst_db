"""
An internal module to download the ESA burst database
"""

import os
import shutil
import subprocess
import tempfile
import zipfile

ESA_DB_URL = "https://sar-mpc.eu/files/S1_burstid_20220530.zip"


def get_esa_burst_db(output_path="burst_map_IW_000001_375887.sqlite3"):
    """Download the ESA burst database."""
    print(f"Downloading ESA burst database from {ESA_DB_URL} to {output_path}.")
    db_filename = "S1_burstid_20220530/IW/sqlite/burst_map_IW_000001_375887.sqlite3"
    cur_dir = os.getcwd()
    output_path = os.path.abspath(output_path)
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            os.chdir(tmpdir)
            subprocess.check_call(["wget", ESA_DB_URL])

            with zipfile.ZipFile(ESA_DB_URL.split("/")[-1], "r") as zip_ref:
                zip_ref.extract(db_filename)
                shutil.move(db_filename, output_path)
                shutil.rmtree(db_filename.split("/")[0])
        finally:
            os.chdir(cur_dir)
