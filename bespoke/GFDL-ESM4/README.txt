GFDL-ESM4
---------

As of 25 May 2023, GFDL-ESM4 has some of its pre-industrial
control (source_id "piControl") and historical (source_id
"historical") simulation data archived in the ESM equivalent
source_ids "esm-piControl" and "esm-hist" (respectively) [1]:

tas for piControl stored under esm-piControl
hfbasin for historical stored under esm-hist

Since very few models (including explicitly Earth-system models)
archive their data under these variants, and there is nothing in
the model documentation paper [2] that suggests separate runs
have been carried out with and without interactive carbon, I
have assumed that these are actually coming from the same
simulation and can thus be mixed.

The scripts in this directory just overwrite the experiment_id
metadata in outputfiles for tas as esm-piControl and hfbasin
as esm-hist (as above) to record this information, but the files
are saved under the same directory/filename structure (e.g.,
tas_yr_piControl_GFDL-ESM4.nc) to simplify later analysis.

References
----------

[1] Eyring, V. et al., 2016: Overview of the Coupled Model
    Intercomparison Project Phase 6 (CMIP6) experimental design
    and organisation, Geosci. Model Dev., 9, 1937-1958,
    doi:10.5194/gmd-9-1937-2016

[2] Dunne, J. P. et al., 2020: The GFDL Earth System Model
    Version 4.1 (GFDL-ESM 4.1): Overall coupled model
    description and simulation characteristics, J. Adv. Model.
    Earth Syst., 12, e2019MS002015, doi:10.1029/2019MS002015
