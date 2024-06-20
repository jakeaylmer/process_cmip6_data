This directory contains simple text files that are read by
src/metadata.py and some bash scripts to determine paths to
various directories, e.g., CMIP6 raw data and processed/output
data. Enter single lines of text without any shell variables,
e.g., "/home/users/myusername/mydirectory/CMIP6",
NOT "~/mydirectory/CMIP6". The file names are assumed by
metadata.py. Ensure the line ending of these files is set to
UNIX (LF), NOT Windows (CR LF), in order for the bash scripts
(e.g., sea_ice/calc_iel.sh) to work.
