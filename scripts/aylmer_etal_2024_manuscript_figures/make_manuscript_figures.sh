#!/bin/bash

# Bash script to generate one or more Aylmer et al. (2024; [1])
# manuscript figures. Figures should be generated from here
# rather than the Python scripts directly, as various settings
# are saved here via hardcoded input arguments to those scripts.
# 
# This script needs to be run from the directory containing the
# script, the various python scripts, and the plotting "utils"
# directory.
# 
# 
# Usage examples
# --------------
# To make and save all figures (14 total):
# 
#     $ bash ./make_manuscript_figures.sh all
# 
# To make and save figures 1, 2, S3, and S6:
# 
#     $ bash ./make_manuscript_figures.sh save 1 2 S3 S6
# 
# To make and show (without saving) figure 3:
# 
#     $ bash ./make_manuscript_figures.sh show 3
# 
# 
# References
# ----------
# [1] Aylmer, J. R., D. Ferreira, and D. Feltham, 2024: Impact
#     of ocean heat transport on sea ice captured by a simple
#     energy balance model, Commun. Earth Environ., accepted in
#     principle May 2024.
# 
# ============================================================ #

if [[ "$*" == *"save"* || "$*" == *"all"* ]]
then
    pyOpts="-u"
    pySaveFig="--savefig"
else
    pyOpts="-i"
    pySaveFig=""
fi

if [[ "$*" == *"all"* ]]
then
    cmdargs=" 1 2 3 4 5 6 7 8 S1 S2 S3 S4 S5 S6"
else
    cmdargs=" $* "
fi

if [[ "${cmdargs}" == *" 1"* ]]
then
    python ${pyOpts} ./ice_edge_vs_oht_and_tas_historical.py \
        --a-xlim 0 4 --a-xticks 0 4 1 \
        --a-ylim -1 5 --a-yticks 0 4 2 \
        --b-xlim -50 110 --b-xticks -50 100 50 \
        --b-ylim -1 5 --b-yticks 0 4 2 \
        --a-cbar-cats -60 60 20 --a-cbar-ticks -60 60 20 \
        --b-cbar-cats 0 3.5 0.5 --b-cbar-ticks 0 3 1 \
        --a-dr-ebm-text 0.0 -0.9 \
        --a-dr-cmip6-text 0.3 -0.5 \
        --a-dr-obs-text -1.1 0.6 \
        --b-dr-ebm-text 0.0 -0.9 \
        --b-dr-cmip6-text 12.0 -1.6 \
        --b-dr-obs-text nan nan \
        --savefigname fig1 ${pySaveFig} \
        --savefigtitle "Historical Arctic sea ice loss, ocean heat transport, and temperature change in CMIP6 and observations"
fi

if [[ "${cmdargs}" == *" 2"* ]]
then
    python ${pyOpts} ./ebm_schematic_maps.py \
        --reflats 65.0 90.0 --savefigname fig2 ${pySaveFig} \
        --savefigtitle "Energy-balance model schematic for the Arctic"
fi

# Common options for Fig. 3 and Fig. S6:
pyScriptOpts="--reflats -60.0 -90.0"
pyScriptOpts="${pyScriptOpts} --a-xlim -0.5 1.5"
pyScriptOpts="${pyScriptOpts} --a-xticks -0.5 1.5 0.5"
pyScriptOpts="${pyScriptOpts} --a-ylim -1 2.5"
pyScriptOpts="${pyScriptOpts} --a-yticks -1 2 1"
pyScriptOpts="${pyScriptOpts} --b-xlim -50 100 "
pyScriptOpts="${pyScriptOpts} --b-xticks -50 100 50"
pyScriptOpts="${pyScriptOpts} --b-ylim -1 2.5"
pyScriptOpts="${pyScriptOpts} --b-yticks -1 2 1"
pyScriptOpts="${pyScriptOpts} --a-cbar-cats -60 60 20"
pyScriptOpts="${pyScriptOpts} --a-cbar-ticks -60 60 20"
pyScriptOpts="${pyScriptOpts} --b-cbar-cats 0 1.5 0.3"
pyScriptOpts="${pyScriptOpts} --b-cbar-ticks 0 1.5 0.3"
pyScriptOpts="${pyScriptOpts} --a-dr-ebm-text -0.25 0.3"
pyScriptOpts="${pyScriptOpts} --a-dr-cmip6-text 0.25 -0.3"
pyScriptOpts="${pyScriptOpts} --a-dr-obs-text -0.35 0.3"
pyScriptOpts="${pyScriptOpts} --b-dr-ebm-text 0.0 0.42"
pyScriptOpts="${pyScriptOpts} --b-dr-cmip6-text 0.0 0.7"
pyScriptOpts="${pyScriptOpts} --b-dr-obs-text 0.0 -0.55"

if [[ "${cmdargs}" == *" 3"* ]]
then
    python ${pyOpts} ./ice_edge_vs_oht_and_tas_historical.py \
        --savefigname fig3 ${pyScriptOpts} ${pySaveFig} \
        --savefigtitle "Historical Southern Ocean sea ice loss, ocean heat transport, and temperature change in CMIP6 and observations"
fi

if [[ "${cmdargs}" == *" S6"* ]]
then
    python ${pyOpts} ./ice_edge_vs_oht_and_tas_historical.py \
        --savefigname figS6 ${pyScriptOpts} ${pySaveFig} \
        --exclude-from-ebm CNRM-CM6-1-HR MPI-ESM1-2-LR \
        --savefigtitle "Historical Southern Ocean sea ice loss, ocean heat transport, and temperature change in sub-sampled CMIP6 and observations"
fi

if [[ "${cmdargs}" == *" 4"* ]]
then
    python ${pyOpts} ./ice_edge_vs_oht_future.py \
        --reflatsn 65.0 90.0 --reflatss -60.0 -90.0 \
        --yravg1 1980 2000 --yravg2 2030 2050 \
        --a-xlim -80 160 --a-xticks -80 160 80 \
        --a-ylim -1 12 --a-yticks 0 12 4 \
        --b-xlim -90 120 --b-xticks -60 120 60 \
        --b-ylim -0.1 2.6 --b-yticks 0 2.5 0.5 \
        --a-dr-ebm-text 0.0 -2.0 \
        --b-dr-ebm-text 0.0 -0.6 \
        --savefigname fig4 ${pySaveFig} \
        --savefigtitle "Future sea ice loss and ocean heat transport change in CMIP6 SSP3-7.0 projections"
fi

if [[ "${cmdargs}" == *" 5"* ]]
then
    python ${pyOpts} ./ebm_parameters_from_piControl.py \
        --a-xlim -2 1.2 --a-xticks -2 2 1 \
        --a-ylim -10 8 --a-yticks -10 10 5 \
        --b-xlim -0.8 1.2 --b-xticks -0.6 1.2 0.6 \
        --b-ylim -6 10 --b-yticks -5 10 5 \
        --c-xlim -8 5 --c-xticks -8 8 4 \
        --c-ylim -2.2 2 --c-yticks -3 2 1 \
        --d-xlim -3 4 --d-xticks -4 4 2 \
        --d-ylim -1 2 --d-yticks -1 2 1 \
        --e-xlim -2.2 2 --e-xticks -3 2 1 \
        --e-ylim -2 1.6 --e-yticks -2 2 1 \
        --f-xlim -0.8 1.2 --f-xticks -0.8 0.8 0.8 \
        --f-ylim -2 4 --f-yticks -2 4 2 \
        --savefigname fig5 ${pySaveFig} \
        --savefigtitle "Parameters for energy-balance model from CMIP6 pre-industrial control simulations"
fi

if [[ "${cmdargs}" == *" 6"* ]]
then
    python ${pyOpts} ./ebm_parameters_variable.py \
        --yravg1 1980 2000 --yravg2 2001 2021 \
        --a-xlim -50 110 --a-xticks -40 80 40 \
        --a-ylim -80 40 --a-yticks -80 40 40 \
        --b-xlim -50 90 --b-xticks -50 50 50 \
        --b-ylim -100 75 --b-yticks -100 50 50 \
        --c-xlim -50 110 --c-xticks -40 80 40 \
        --c-ylim -25 40 --c-yticks -20 40 20 \
        --d-xlim -50 90 --d-xticks -50 50 50 \
        --d-ylim -40 80 --d-yticks -40 80 40 \
        --savefigname fig6 ${pySaveFig} \
        --savefigtitle "Parameters for energy-balance model from CMIP6 historical simulations"
fi

if [[ "${cmdargs}" == *" 7"* ]]
then
    python ${pyOpts} ./ebm_parameters_variable.py \
        --yravg1 1980 2000 --yravg2 2030 2050 \
        --a-xlim -70 140 --a-xticks -40 120 40 \
        --a-ylim -110 100 --a-yticks -100 100 50 \
        --b-xlim -90 120 --b-xticks -60 120 60 \
        --b-ylim -150 100 --b-yticks -150 100 50 \
        --c-xlim -70 140 --c-xticks -40 120 40 \
        --c-ylim -30 60 --c-yticks -30 60 30 \
        --d-xlim -90 120 --d-xticks -60 120 60 \
        --d-ylim -50 150 --d-yticks -50 150 50 \
        --savefigname fig7 ${pySaveFig} \
        --savefigtitle "Parameters for energy-balance model from CMIP6 SSP3-7.0 simulations"
fi

if [[ "${cmdargs}" == *" 8"* ]]
then
    python ${pyOpts} ./oht_observation_estimate_method.py \
        --a-xlim 1978 2023 --a-xticks 1980 2020 10 \
        --a-ylim 225 280 --a-yticks 220 280 20 \
        --b-xlim 1978 2023 --b-xticks 1980 2020 10 \
        --b-ylim 290 430 --b-yticks 300 420 40 \
        --savefigname fig8 ${pySaveFig} \
        --savefigtitle "Estimating the Circulation and Climate of the Ocean (ECCO) 1980 to 2021 ocean heat transport change estimation"
fi

if [[ "${cmdargs}" == *" S1"* ]]
then
    python ${pyOpts} ./ice_edge_vs_gmst_historical.py \
        --a-xlim 0.15 1.1 --a-xticks -0.2 1 0.2 \
        --a-ylim -1 5 --a-yticks 0 4 2 \
        --b-xlim 0.15 1.1 --b-xticks -0.2 1 0.2 \
        --b-ylim -0.4 1.7 --b-yticks -0.4 1.6 0.4 \
        --a-dr-cmip6-text -0.05 -1.6 \
        --a-dr-obs-text -0.2 0.55 \
        --b-dr-cmip6-text -0.4 -0.1 \
        --b-dr-obs-text nan nan \
        --savefigname figS1 ${pySaveFig} \
        --savefigtitle "Historical sea ice loss and global mean temperature change in CMIP6 and observations"
fi

if [[ "${cmdargs}" == *" S2"* ]]
then
    python ${pyOpts} ./lagged_correlations_ensemble.py \
        --a-xlim -25 25 --a-xticks -20 20 10 \
        --savefigname figS2 ${pySaveFig} \
        --savefigtitle "Lagged correlations between polar variables in the CMIP6 ensemble"
fi

if [[ "${cmdargs}" == *" S3"* ]]
then
    python ${pyOpts} ./lagged_correlations_individual_models.py \
        --a-xlim -25 25 --a-xticks -20 20 10 \
        --savefigname figS3 ${pySaveFig} \
        --savefigtitle "Lagged correlations between polar variables in selected CMIP6 models"
fi

if [[ "${cmdargs}" == *" S4"* ]]
then
    python ${pyOpts} ./ebm_schematic_maps.py \
        --reflats -60.0 -72.0 --savefigname figS4 ${pySaveFig} \
        --savefigtitle "Energy-balance model schematic for the Southern Ocean"
fi

if [[ "${cmdargs}" == *" S5"* ]]
then
    python ${pyOpts} ./ice_edge_vs_oht_and_tas_historical.py \
        --reflats -60.0 -72.0 \
        --a-xlim -1 2 --a-xticks -1 2 1 \
        --a-ylim -1 2.5 --a-yticks -1 2 1 \
        --b-xlim -2 4 --b-xticks -2 4 2 \
        --b-ylim -1 2.5 --b-yticks -1 2 1 \
        --a-cbar-cats -3 3 1 --a-cbar-ticks -3 3 1 \
        --b-cbar-cats 0 1.5 0.3 --b-cbar-ticks 0 1.5 0.3 \
        --a-dr-ebm-text -0.35 0.25 \
        --a-dr-cmip6-text 0.25 -0.3 \
        --a-dr-obs-text -0.5 0.3 \
        --b-dr-ebm-text 0.0 0.4 \
        --b-dr-cmip6-text 0.0 0.7 \
        --b-dr-obs-text 0.0 -0.55 \
        --savefigname figS5 ${pySaveFig} \
        --savefigtitle "Historical Southern Ocean sea ice loss, ocean heat transport convergence, and temperature change in CMIP6 and observations"
fi

exit 0
