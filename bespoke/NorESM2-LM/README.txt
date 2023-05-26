NorESM2-LM / NorESM2-MM
-----------------------

For these models, the ocean grid has one extra row of data on
axis = 1 (y direction) compared to the sea ice data. The model
documentation explains the issue and that for sea ice data the
final row of any corresponding ocean grid variable (here
areacello) should be dropped.

All other ocean and atmospheric diagnostics can be processed
using the general routines.

Reference:

https://noresm-docs.readthedocs.io/en/noresm2/faq/postp_plotting_faq.html

[Section 4.2 "Different sea-ice and ocean grid", accessed 24/05/2023]
