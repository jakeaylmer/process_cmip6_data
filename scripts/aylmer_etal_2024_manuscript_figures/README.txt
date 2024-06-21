This directory contains scripts that generate the eight main
figures and six supplementary figures in the work by Aylmer et
al. (2024) [1]. In the current directory (of this README.txt),
there are nine Python scripts:

ebm_parameters_from_piControl.py         --> Fig. 5
ebm_parameters_variable.py               --> Figs. 6, 7
ebm_schematic_maps.py                    --> Figs. 2, S4
ice_edge_vs_gmst_historical.py           --> Fig. S1
ice_edge_vs_oht_and_tas_historical.py    --> Figs. 1, 3, S5, S6
ice_edge_vs_oht_future.py                --> Fig. 4
lagged_correlations_ensemble.py          --> Fig. S2
lagged_correlations_individual_models.py --> Fig. S3
oht_observation_estimate_method.py       --> Fig. 8

However, it is easier (and intended) to use the bash script,
make_manuscript_figures.sh, as this contains the various options
to each Python script to generate each figure as it appears in
the manuscript, and parses command-line arguments to show or
save specific figure(s). See docs of that script, but, e.g.:

    $ bash ./make_manuscript_figures.sh show 5

will show Fig. 5 interactively, while:

    $ bash ./make_manuscript_figures.sh save 3 S5 S6

will save Figs. 3, S5, and S6 to file (.svg and .png formats).
Only one figure can be "show"n at a time, but multiple can be
"save"d in one use of the script. Note that it should be run
from this directory.

The scripts above are predominantly just calls to the multitude
of shared routines in the subdirectory "utils":

    utils/ebm.py            calculation of EBM parameters/slopes
    utils/maths.py          general mathematical routines
    utils/model_lists.py    which models to load per experiment
    utils/observations.py   loads/prepares observations
    utils/plot_style.py     plot appearance, fonts/element sizes
    utils/plot_tools.py     common plot routines, figure layouts
    utils/script_tools.py   loading cmip6 processed data,
                            command-line argument parsing


Plot style: fonts
-----------------
Figures use "Nimbus Sans" fonts, which is available on the
author's local system and visible to matplotlib automatically.
If this font is not available, matplotlib will select the next
available font (most likely, "DejaVu Sans", the default font
packaged with matplotlib) and plots will thus look a little
different (and possibly the font size, utils.plot_style.fs_0,
would need adjusting to accommodate a different font so that
labels do not overlap with other plot elements).



Paths
-----
The utils code uses the paths directory (i.e., that used by
the processing scripts) to locate data and save figures. The
following paths files should be set for these scripts to work:

    paths/path_atmospheric_reanalyses_processed_data.txt
    paths/path_cmip6_processed_data.txt
    paths/path_nsidc_processed_data.txt

which are shared anyway by the processing code. Note that the
data loading functions in utils/script_tools.py using the "api"
of the processing code, so this must also be visible on the
Python PATH.

Also, two additional paths are set in that directory in the same
format (i.e., single line giving full path to file or
directory). Firstly:

    paths/path_ecco_oht.txt

should point to the ECCO OHT data obtained from
doi:10.5281/zenodo.7869067; specifically, it should point to the
file "MHT.jld2" contained therein (not just the directory that
contains it; this is unlike the other paths which just point to
top-level directories because this is just data that we have
downloaded and not applied any additional processing to).

Secondly, the path set in:

    paths/path_manuscript_figures.txt

determines where saved figures go. The directory is created if
it does not exist, and leaving it blank sets it to just save
figures to the current working directory.



Library versions
----------------
The version of Python used/tested with is 3.12.0 on Linux. The
modules required, and the versions tested with, are:

    cartopy      0.22.0   (only in ebm_schematic_maps.py)
    numpy        1.26.0
    netCDF4      1.6.5
    matplotlib   3.8.1
    scipy        1.11.3
    tabulate     0.9.0



Further processing of figure files
----------------------------------
There is another bash script, convert_svg_to_pdf.sh, which was
just used to produce final versions submitted to the publisher
in PDF format (there is a note in that script why direct PDF
output from matplotlib was not used).

There is a special case with Fig. S4 (Southern Ocean EBM
schematic): the figure is produced upside down. This is due to a
limitation with the Python module cartopy (or, indeed, the
author's expertise with it). The output was rotated through 180
degrees manually using Inkscape before exporting to the final
PDF version.



References
----------
[1] Aylmer, J. R., D. Ferreira, and D. L. Feltham, 2024: Impact
    of ocean heat transport on sea ice captured by a simple
    energy balance model, Commun. Earth Environ., accepted in
    principle May 2024.
