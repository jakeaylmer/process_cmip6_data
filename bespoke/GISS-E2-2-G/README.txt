GISS-E2-2-G
-----------

The data for hfds for this model is provided on the atmosphere
grid -- presumed, because the data dimensions match and also the
siconc data is only available on the atmosphere grid (siconca).
There is no obvious documentation verifying this, and the CMIP6
attribute "grid_label" is set to "gn" (which implies native,
ocean grid) rather than "gr" (which is what the data implies). A
separate script is provided for saving hfds yearly fields which
loads the atmosphere grid data rather than the ocean for saving.

Regardless of realm, the ocean grid for this model has
independent longitude and latitude axes, which means that the
area integrals/averages can be computed exactly from the
coordinate bounds (like the atmosphere diagnostics), so a
separate script is provided for this too.

