AWI-CM-1-1-MR
-------------

This model uses an unstructured ocean grid, and ocean variables
are provided as 1D arrays (i.e., a single list of all grid
cells). These scripts accommodate this unique structure (for
hfds and opottemptend and associated diagnostics, and sea ice
area/extent).

The ocean grid cell area (areacello) is not saved in standard
format for this model (as it cannot be) and the above scripts
load the raw areacello data.

For the sea ice-edge latitude, data does not have any land-mask
information and does not respond well to bilinear interpolation.
Scripts are provided that interpolates data to a regular 0.25
degree rectilinear grid using distance-weighted remapping, and
a land mask is added to the remapped data.

Note also that sea ice concentration data is only available
daily (these should be converted to monthly-mean fields before
processing).

The atmospheric diagnostics can be processed using the general
routines.
