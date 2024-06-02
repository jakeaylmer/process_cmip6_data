#!/bin/bash

# This script runs directly on the raw netCDF files for the
# passive microwave sea ice concentration data downloaded from
# the National Snow and Ice Data Center (NSIDC). Specifically,
# monthly data from the latest versions of the two SMMR and
# SSM/I-SSMIS datasets [1,2] commonly referred to as "NASA Team"
# and "Bootstrap" sea ice concentrations.
# 
# The main purpose for including here is to make available the
# process for obtaining observational estimates of the change in
# Arctic and Antarctic sea ice edge over the satellite era used
# in the work by Aylmer et al. 2024 [3], for which this code
# repository is primarily for. While the Python code in this
# directory does use some of the src modules for convenience the 
# procedure is mostly independent of the CMIP6/atmospheric
# reanalysis routines. This is because additional processing
# steps are required, particularly in the northern hemisphere
# which requires masking of the satellite pole hole (which
# changes size with changes in instrument/satellite source)
# and application of a 'valid ice' mask. The ice edge latitude
# calculation is similar to that for the CMIP6 models, but that
# uses external routines anyway [4].
# 
# This bash script carries out all of the required steps for
# both datasets and all years, each of which uses a separate
# Python script except for step (3), the interpolation, which
# uses CDO and is done in-line here: 
# 
#     (1) Prepare raw data into a format that can be
#         interpolated using CDO. This intermediate data is
#         saved to a temporary directory and deleted at the end.
#     (2) Saves the climatology sea ice concentration data used
#         in the maps (on the native grid; Fig. 2 and
#         supplementary Fig. 4 of [3])
#     (3) Interpolates the concentration data to a regular 0.5
#         degree lon/lat grid (not global as there is no need
#         and anyway the raw data is not global; see settings
#         below)
#     (4) Creates pole hole mask auxilliary data for the north-
#         ern hemisphere (on both native and interpolated grids;
#         these masks are required for steps 5 and 6)
#     (5) Calculates and saves sea ice area and extent
#         diagnostics, on both native and interpolated grids for
#         area/extent. Note these diagnostics are not explicitly
#         used in [3].
#     (6) Calculates sea ice-edge latitude diagnostics (monthly
#         and annual means, zonal and zonal-means).
# 
# References
# ----------
# [1] DiGirolamo, N., C. L. Parkinson, D. J. Cavalieri, P.
#     Gloersen, and H. J. Zwally, 2022: Sea Ice Concentrations
#     from Nimbus-7 SMMR and DMSP SSM/I-SSMIS Passive Microwave
#     Data, Version 2 [Data Set], Boulder, Colorado USA, NASA
#     National Snow and Ice Data Center Distributed Active
#     Archive Center, doi:10.5067/MPYG15WAA4WX
# 
# [2] Comiso, J. C., 2023: Bootstrap Sea Ice Concentrations from
#     Nimbus-7 SMMR and DMSP SSM/I-SSMIS, Version 4 [Data Set],
#     Boulder, Colorado USA, NASA National Snow and Ice Data
#     Center Distributed Active Archive Center,
#     doi:10.5067/X5LG68MH013O
# 
# [3] Aylmer, J. R., D. Ferreira, and D. L. Feltham, 2024:
#     Impact of ocean heat transport on sea ice captured by a
#     simple energy balance model, Commun. Earth Environ.,
#     accepted in principle May 2024
# 
# [4] Aylmer, J. R., 2021: Sea ice-edge latitude diagnostic
#     code, doi:10.5281/zenodo.5494523
# ============================================================ #


# ------------------------------------------------------------ #
# Settings/paths
# ============================================================ #

# Flags to switch on/off each step (obviously, intermediate data
# from prior steps must still exist; just allows re-running of
# later steps). Set to strings "true" or "false":
doStepPreparation="true"  # step (1)
doStepClimatology="false"  # step (2)
doStepInterpolate="true"  # step (3)
doStepMakeNPmasks="false"  # step (4)
doStepAreaExtents="false"  # step (5)
doStepIceEdgeLats="true"  # step (6)

keepIntermediateData="true"

# Start and end years of data, inclusive (whole years assumed):
ys=1979
ye=2022

# Target grid definitions for interpolation as interpreted as a
# text file by CDO genbil, one each for northern and southern
# hemispheres (the following written to a temporary text file)
# These are regular 0.5 degree grids (not global):
tgridn="gridtype = lonlat\n"
tgridn="${tgridn}xsize    = 720\n"
tgridn="${tgridn}ysize    = 120\n"
tgridn="${tgridn}xfirst   = -179.75\n"
tgridn="${tgridn}xinc     = 0.5\n"
tgridn="${tgridn}yfirst   = 30.25\n"
tgridn="${tgridn}yinc     = 0.5\n"

tgrids="gridtype = lonlat\n"
tgrids="${tgrids}xsize    = 720\n"
tgrids="${tgrids}ysize    = 110\n"
tgrids="${tgrids}xfirst   = -179.75\n"
tgrids="${tgrids}xinc     = 0.5\n"
tgrids="${tgrids}yfirst   = -89.75\n"
tgrids="${tgrids}yinc     = 0.5\n"

# Which years to use for generating each pole mask (see time
# validity ranges in script make_pole_hole_mask.py - otherwise,
# it doesn't matter, it's just to see where data is missing
# around the north pole for the different instrument periods):
yearPoleMaskSMMR=1985
yearPoleMaskSSMI=2005
yearPoleMaskSSMIS=2015

# Locations of raw data (each contains subdirectories YYYY
# with the files for year YYYY exactly as downloaded from the
# NSIDC):
_rawDataDir="/storage/basic/cpom/gb919150/NSIDC"
rawData0051N="${_rawDataDir}/NSIDC-0051_nasateam_v2/raw/nh/monthly/"
rawData0051S="${_rawDataDir}/NSIDC-0051_nasateam_v2/raw/sh/monthly/"
rawData0079N="${_rawDataDir}/NSIDC-0079_bootstrap_v4/raw/nh/monthly/"
rawData0079S="${_rawDataDir}/NSIDC-0079_bootstrap_v4/raw/sh/monthly/"

# Directory containing this bash script and related Python
# scripts:
scriptsDir="${HOME}/phd/process_cmip6_data/bespoke/passive_microwave"

# Set a directory for intermediate files (for 1979-2022 monthly
# data for both hemispheres and datasets, uses about 840 MB):
wrkDir="${HOME}/tmpPassiveMicrowaveProcessing"

# ------------------------------------------------------------ #


mkdir -p ${wrkDir}
cd ${wrkDir}

# (1) Prepare sea ice concentration using a python script to
#     combine all monthly files into yearly files, set the
#     missing flags, add the coordinate data, and reverse the
#     axis=1 dimension which is wrong in the raw data.
if [[ "${doStepPreparation}" == "true" ]]
then
    for y in $(seq ${ys} ${ye})
    do
        python -u ${scriptsDir}/prepare_raw_siconc_data.py \
            -i ${rawData0051N}/${y}/*.nc --hemisphere "n" \
            -d "NSIDC-0051" -o ${wrkDir}/NSIDC-0051_n_${y}.nc
        
        python -u ${scriptsDir}/prepare_raw_siconc_data.py \
            -i ${rawData0051S}/${y}/*.nc --hemisphere "s" \
            -d "NSIDC-0051" -o ${wrkDir}/NSIDC-0051_s_${y}.nc
        
        python -u ${scriptsDir}/prepare_raw_siconc_data.py \
            -i ${rawData0079N}/${y}/*.nc --hemisphere "n" \
            -d "NSIDC-0079" -o ${wrkDir}/NSIDC-0079_n_${y}.nc
        
        python -u ${scriptsDir}/prepare_raw_siconc_data.py \
            -i ${rawData0079S}/${y}/*.nc --hemisphere "s" \
            -d "NSIDC-0079" -o ${wrkDir}/NSIDC-0079_s_${y}.nc
    done
fi


# (2) Save climatologies from the above prepared data (in [3]
#     just Bootstrap, 0079, is used in Figs. 2 and S4, but save
#     both anyway):
if [[ "${doStepClimatology}" == "true" ]]
then
    for ds in "NSIDC-0051" "NSIDC-0079"
    do
        for hemi in "n" "s"
        do
            python -u -W ignore \
                ${scriptsDir}/calc_siconc_climatology.py \
                -i ${wrkDir}/${ds}_${hemi}_{1980..2000}.nc \
                --hemisphere ${hemi} -d ${ds}
            
            python -u -W ignore \
                ${scriptsDir}/calc_siconc_climatology.py \
                -i ${wrkDir}/${ds}_${hemi}_{2001..2021}.nc \
                --hemisphere ${hemi} -d ${ds}
        done
    done
fi



# (3) Interpolation:
if [[ "${doStepInterpolate}" == "true" ]]
then
    printf "${tgridn}" > ${wrkDir}/_tgridn.txt
    printf "${tgrids}" > ${wrkDir}/_tgrids.txt
    
    for hemi in "n" "s"
    do
        cdo genbil,${wrkDir}/_tgrid${hemi}.txt \
            ${wrkDir}/NSIDC-0051_${hemi}_2000.nc \
            ${wrkDir}/_weights_${hemi}.nc
        
        for ds in "NSIDC-0051" "NSIDC-0079"
        do
            for y in $(seq ${ys} ${ye})
            do
                cdo remap,${wrkDir}/_tgrid${hemi}.txt,${wrkDir}/_weights_${hemi}.nc \
                    ${wrkDir}/${ds}_${hemi}_${y}.nc \
                    ${wrkDir}/${ds}_${hemi}_${y}_remapped.nc
            done
        done
    done
fi


# (4) Generate pole hole mask auxilliary data:
if [[ "${doStepMakeNPmasks}" == "true" ]]
then
    python -u -W ignore ${scriptsDir}/make_pole_hole_mask.py \
        -i ${wrkDir}/NSIDC-0051_n_${yearPoleMaskSMMR}.nc \
        --whichmask "SMMR" -d "NSIDC-0051"
    
    python -u -W ignore ${scriptsDir}/make_pole_hole_mask.py \
        -i ${wrkDir}/NSIDC-0051_n_${yearPoleMaskSSMI}.nc \
        --whichmask "SSMI" -d "NSIDC-0051"
    
    python -u -W ignore ${scriptsDir}/make_pole_hole_mask.py \
        -i ${wrkDir}/NSIDC-0051_n_${yearPoleMaskSSMIS}.nc \
        --whichmask "SSMIS" -d "NSIDC-0051"
        
    python -u -W ignore ${scriptsDir}/make_pole_hole_mask.py \
        -i ${wrkDir}/NSIDC-0051_n_${yearPoleMaskSMMR}_remapped.nc \
        --whichmask "SMMR" -d "NSIDC-0051"
    
    python -u -W ignore ${scriptsDir}/make_pole_hole_mask.py \
        -i ${wrkDir}/NSIDC-0051_n_${yearPoleMaskSSMI}_remapped.nc \
        --whichmask "SSMI" -d "NSIDC-0051"
    
    python -u -W ignore ${scriptsDir}/make_pole_hole_mask.py \
        -i ${wrkDir}/NSIDC-0051_n_${yearPoleMaskSSMIS}_remapped.nc \
        --whichmask "SSMIS" -d "NSIDC-0051"
fi


# (5) Calculate sea ice area and extent:
if [[ "${doStepAreaExtents}" == "true" ]]
then
    
    for hemi in "n" "s"
    do
        for ds in "NSIDC-0051" "NSIDC-0079"
        do
            python -u -W ignore \
                ${scriptsDir}/calc_sia_sie.py \
                -i ${wrkDir}/${ds}_${hemi}_????.nc \
                -d ${ds} --hemisphere ${hemi}
            
            python -u -W ignore \
                ${scriptsDir}/calc_sia_sie.py \
                -i ${wrkDir}/${ds}_${hemi}_????_remapped.nc \
                -d ${ds} --hemisphere ${hemi}
        done
    done
fi


# (6) Calculate sea ice-edge latitudes:
if [[ "${doStepIceEdgeLats}" == "true" ]]
then
    
    for hemi in "n" "s"
    do
        for ds in "NSIDC-0051" "NSIDC-0079"
        do
            python -u -W ignore \
                ${scriptsDir}/calc_iel.py \
                -i ${wrkDir}/${ds}_${hemi}_????_remapped.nc \
                -d ${ds} --hemisphere ${hemi}
        done
    done
fi

if [[ "${keepIntermediateData}" != "true" ]]
then
    echo "rm -dr ${wrkDir}"
    rm -dr ${wrkDir}
fi

exit 0
