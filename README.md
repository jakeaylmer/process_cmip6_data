[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12207486.svg)](https://doi.org/10.5281/zenodo.12207486)

# process_cmip6_data

This repository contains the code used to process data from the CMIP6 archive into the diagnostics required in the work by Aylmer et al. (2024).[^manuscript] It is structured as a Python package with numerous Python (and some bash) scripts included in various subdirectories. It may be useful for/adaptable to processing other CMIP diagnostics. Processed data is saved in NetCDF format, with CF-compliant (as far as possible) metadata including full citations to the source data and relevant model descriptions. All of the code and scripts are documented individually with their specific usage; this README just gives a 'high-level' overview of the steps used to generate the data appearing in the aforementioned study. All processed data are available online.[^data] Scripts are also included in this repository that reproduce the manuscript[^manuscript] figures.


## Usage

### Repository structure overview

There are two 'code' directories (sub-packages):
  * `/src`: the main data processing code
  * `/api`: provides an interface for loading data and is imported when the package is imported

Then there is a 'configuration' directory, `/paths`, containing plain-text files, each of which is a single line pointing to a specific directory (e.g., raw data; see _Raw data and paths_). The remaining directories contain scripts:
  * `/atmosphere`, `/ocean`, and `/sea_ice`: main/generic processing scripts for computing diagnostics in the atmosphere, ocean, and sea ice domains, respectively
  * `/bespoke/*`: modified copies of the main processing scripts for specific cases where there is a slight difference in the raw data structure or some other issue that needs a 'bespoke' approach (see _'Bespoke' processing scripts_)
  * `/qc`: 'quality-control' plots (see _'Quality control' scripts_)
  * `/scripts`: other/miscellaneous scripts (see _Other scripts_)


### Raw data and paths
Raw data is not provided in the repository but is publicly available from the CMIP6 archive, accessible from a number of data hosts (see, e.g., [https://pcmdi.llnl.gov/CMIP6/](https://pcmdi.llnl.gov/CMIP6/)). This data should first be downloaded into some directory, which can be called anything, but for example `/path/to/raw/cmip6/data`. This 'top-level' directory should have a subdirectory structure `./source_id/experiment_id/variable_id` (if [downloading CMIP6 data using](https://esgf.github.io/esgf-user-support/user_guide.html) `wget`, this can be set using `&download_structure=source_id,experiment_id,variable`). An example would be `./UKESM1-0-LL/piControl/siconc/` containing all NetCDF files for that particular model, experiment, and variable (in this case, sea ice concentration for the pre-industrial control simulation of the model UKESM1-0-LL). Multiple ensemble members are included in the same directory.

File names of raw data should not be modified; there is a standard naming convention used by the CMIP6 data archive that is assumed by the `/src/load_raw_data.py` module. However, that module does need to know where this data is: in the directory `/paths`, there is a plain text file `path_cmip6_raw_data.txt`, which needs to contain the 'top-level' directory&#x2014;in the example above, "/path/to/raw/cmip6/data"&#x2014;before any processing scripts are run. The scripts accept command-line arguments to specify the model (`-m UKESM1-0-LL`) and experiment (`-x piControl`), from which the code in `/src/load_raw_data.py` finds and loads the relevant raw data, concatenating in time as required.

Similarly, the file `/paths/path_cmip6_processed_data.txt` needs to be edited with the desired output path for processed diagnostics (the subdirectory structure is determined automatically).

There are a few other text files in `/paths` used for similar&#x2014;and in most cases self-explanatory&#x2014;purposes.

Note that the path of the directory containing the root directory respository needs to be on the Python `PATH` environment variable, so that it is importable as a module within Python scripts (e.g., on Unix systems: `export PYTHONPATH=${PYTHONPATH}:/path/to/directory`).


### Source code
There main processing code is under the `/src` directory. The scripts (see _Script directories_) are mostly just calls to these functions/routines:
  * `/src/diagnostics.py`
      * general diagnostics (e.g., time averaging, spatial averaging/integrals)
  * `/src/load_processed_data.py`
      * for loading data that has been processed already
  * `/src/load_raw_data.py`
      * for loading raw CMIP6 data
  * `/src/metadata.py`
      * some default settings included here, but mainly defines metadata for each model/variable. The 20 models analysed by Aylmer et al. (2024)[^manuscript] are already defined here; other models would need to be added separately before the processing scripts would work.
      * also includes similar kinds of metadata for atmospheric reanalyses and passive microwave data but these are only accessed by the relevant `/bespoke/*` scripts (see _Script directories_)
  * `/src/netcdf.py`
      * writing processed data to netCDF
      * default variable/dimension names, global attributes, and file name formats are also defined here
  * `/src/qc.py`
      * 'Quality-control' plotting code; uses the api loading routines (thus also acting as a test for that; see _API_)
  * `/src/script_tools.py`
      * parsing command-line arguments
  * `/src/utils.py`
      * general Python functions used throughout


### API
The `/api` sub-package is provided to make loading data easier because there are a large number of diagnostic outputs and the file names which, while descriptive, are verbose. In some cases there are multiple methods for the same diagnostics (ocean heat transport being the notorious example) and so to load this data for all models would require keeping track of which method is used for each model. This API implements a set of unique "keywords" and "aliases" for each diagnostic. The keywords are short and simple, for example "oht", and the API routines accept these, specifically via the module `/api/load_data.py`, loading the appropriate diagnostic for each specified model. The "aliases" are shortened versions of the full diagnostic names and are primarily used internally, though can also be specified in case multiple methods are used for one model.

The API also uses the main `/src` code (`metadata.py`, `netcdf.py`, and `utils.py` modules) to determine file names (so if file names are changed after being saved with the processing code, the API will not work) and the location of processed data (which is set by `/paths/path_cmip6_processed_data.txt` and _can_ be changed after processing if the processed data is moved). 

For example usage of the API, see the manuscript figure generation scripts. For a list of all "keywords" and "aliases" that have been defined, the script `/scripts/print_defined_diagnostics.py` is provided.


### Script directories

#### Main CMIP6 processing scripts
These are found under the directories `/atmosphere`, `/ocean`, and `/sea_ice`, and are used to compute several diagnostics in each domain. They apply to any model and experiment, changed via the command line arguments, e.g., `-x historical -m CESM2`. What each script does is hopefully self-explanatory from its file name (e.g., `/atmosphere/calc_tas_area_averages.py`), but is also briefly documented at the beginning of each file. In some cases, one script needs to be run after another which generates some intermediate data (also documented where applicable).


#### 'Quality control' scripts
Under the directory `/qc`, there are a number of scripts that generate quick (i.e., minimally formatted) plots of the processed data.


#### 'Bespoke' processing scripts
Sometimes, CMIP6 data are inconsistent and the main processing scripts need adjustments&#x2014;for example, the model AWI-CM-1-1-MR uses an unstructured ocean grid, requiring special treatment for ocean and sea ice diagnostics. This directory contains copies of those scripts modified slightly, for one reason or another, for specific models and diagnostics. There is a subdirectory for each model that we analysed that required this, each of which includes a `README.txt` explaining what the issues were and what modifications have been made.


##### Atmospheric reanalyses
In our work, we also analysed near-surface temperature in atmospheric reanalyses. Although these are not CMIP6 data, we used the same code and workflow as that for the CMIP6 data and thus include the relevant scripts in the subdirectory `/bespoke/atmospheric_reanalyses`. The various data and model citations are included in `src/metadata.py` and embedded in the netCDF outputs. The processed datasets are available with the archived CMIP6 processed data.[^data]


##### Passive microwave sea ice concentration data
We also required observations of sea ice, so include the scripts here under `/bespoke/passive_microwave`. The raw data comes from passive microwave observations of sea ice concentration obtained from the National Snow and Ice Data Center (NSIDC).[^nsidc0051][^nsidc0079] There are significantly more steps required to process this data (documented in the code) and so these scripts are largely independent of the `/sea_ice` scripts for CMIP6; however, they do make use of the `/src` code, including data citations in `/src/metadata.py` and using some of the functions in `/src/netcdf.py`. The processed datasets are available with the archived CMIP6 processed data.[^data]



#### Other scripts
The `/scripts` directory contains some other scripts:

  * `/scripts/print_defined_diagnostics.py`: this prints the "keywords" and "aliases" defined by the API code
  * `/scripts/update_global_nc_attributes.py`: this ad-hoc script was used to update the metadata of all processed CMIP6 diagnostics upon acceptance of the manuscript (it does not modify the data)
  * `/scripts/aylmer_etal_2024_manuscript_figures`: this directory contains scripts and code generating our manuscript figures. See `README.txt` in that directory for details and usage.


## Software versions

These are the versions of software and libraries used and thus known to work with all code in this repository:
  * Operating system: Linux (Rocky Linux 8.10; [https://rockylinux.org/](https://rockylinux.org/))
      * Linux is required due to the use of Bash scripts and CDO (see below) in the calculation of some diagnostics, notably the sea ice-edge latitude. However, those which are purely in Python are likely to work on other operating systems but this has not been tested.
  * Bash 4.4.20(1)-release
  * Climate Data Operators (CDO) version 2.3.0 ([https://mpimet.mpg.de/cdo](https://mpimet.mpg.de/cdo))
  * Python version 3.12.0 ([https://www.python.org/downloads/release/python-3120/](https://www.python.org/downloads/release/python-3120/)), with the following packages:
      * **Required for processing data:**
          * `netCDF4` version 1.6.5
          * `numpy` version 1.12.0
          * `tabulate` version 0.9.0 (pip installable)
          * `ice_edge_latitude`, version 1.0.0 (available at [doi:10.5281/zenodo.5494524](https://doi.org/10.5281/zenodo.5494524))
      * **Additional requirements for quality control plots and manuscript figures:**
          * `matplotlib` version 3.8.1
      * **Additional requirements for manuscript figures:**
          * `cartopy` version 0.22.0
          * `scipy` version 1.11.3


## Limitations

* Processing code is only setup to handle monthly-frequency input data (except for ocean heat content tendency variables "opottemptend" and "ocontemptend"; the scripts under `/ocean` using these  expects yearly data)
* CMIP6 raw data file names should be as obtained from the CMIP6 archive and not changed
* This repository was developed with our specific applications in mind: no claims are made with regards to the robustness of the Python code, its conformity with PEP 8, its comprehensibility, nor whether the overall structure of the package is optimal


[^manuscript]: Aylmer, J. R., D. Ferreira, and D. L. Feltham, 2024: Impact of ocean heat transport on sea ice captured by a simple energy balance model, _Commun. Earth Environ._, **5**, 406, doi:[10.1038/s43247-024-01565-7](https://doi.org/10.1038/s43247-024-01565-7)
[^data]: Aylmer, J. R., 2024: Diagnostics from CMIP6, atmospheric reanalyses, and passive-microwave observations used to examine the impact of ocean heat transport on Arctic and Antarctic sea ice [Data Set], _University of Reading_, [doi:10.17864/1947.001333](https://doi.org/10.17864/1947.001333)
[^nsidc0051]: DiGirolamo, N., C. L. Parkinson, D. J. Cavalieri, P. Gloersen, and H. J. Zwally, 2022: Sea ice concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS passive microwave data, version 2 [Data Set], Boulder, Colorado USA, NASA National Snow and Ice Data Center Distributed Active Archive Center, [doi:10.5067/MPYG15WAA4WX](https://doi.org/10.5067/MPYG15WAA4WX)
[^nsidc0079]: Comiso, J. C. (2023): Bootstrap sea ice concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS, version 4 [Data Set], Boulder, Colorado USA, NASA National Snow and Ice Data Center Distributed Active Archive Center, [doi:10.5067/X5LG68MH013O](https://doi.org/10.5067/X5LG68MH013O)
