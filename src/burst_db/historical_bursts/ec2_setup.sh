#/usr/bin/env bash
set -e

curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge-$(uname)-$(uname -m).sh -b

source ./mambaforge//etc/profile.d/conda.sh 
source ./mambaforge//etc/profile.d/mamba.sh 
