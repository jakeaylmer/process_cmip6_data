"Bespoke" processing scripts
----------------------------

Directories that are a CMIP6 model name:
    These contain copies of the main processing scripts under
    atmosphere, ocean, and sea_ice, modified ad-hoc for
    specific models that, for some reason or another, need an
    extra step or some slightly different processing. Those
    reasons are documented in a README.txt under each directory.


Directory "atmospheric_reanalyses":
    This directory contains scripts that compute the near-
    surface air temperature diagnostics for atmospheric
    reanalyses (CFSR, CFSv2, ERA5, JRA-55, and MERRA-2). These
    are slightly modified versions of the CMIP6 scripts under
    atmosphere, mainly for slight differences in netCDF
    metadata and different time spans. The src/metadata.py
    module contains a section for reanalyses metadata used by
    these scripts.


Directory "passive_microwave"
    This directory contains scripts that compute sea ice area,
    extent, and sea ice-edge latitude diagnostics from passive-
    microwave observations of sea ice concentration, for two
    datasets provided by the National Snow and Ice Data Center.
    They are largely independent of the CMIP6 sea ice processing
    scripts (i.e., those under sea_ice) but utilise the src/
    netcdf.py and metadata.py modules for convenience. The main
    purpose for including these are to make available the
    procedure used for computing observational estimates of the
    change in sea ice edge over the satellite era in [1];
    otherwise, these don't really "go" with the CMIP6 processing
    code since several additional steps/masks/processing are
    required for the passive microwave data.


References
----------
[1] Aylmer, J. R., D. Ferreira, and D. L. Feltham, 2024: Impact
    of ocean heat transport on sea ice captured by a simple
    energy balance model, Commun. Earth Environ., accepted in
    principle May 2024.
