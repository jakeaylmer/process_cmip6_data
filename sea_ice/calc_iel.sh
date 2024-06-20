#!/bin/bash

# Interpolates raw monthly sea ice concentration data to a
# regular longitude-latitude grid used climate data operators
# (CDO), then running the companion python script "calc_iel.py"
# to calculate sea ice-edge latitude diagnostics from the
# interpolated data.
# 
# Command line options are available, defaults of which are:
model="CESM2"             # model_id
experiment="piControl"    # experiment_id
members=("r1i1p1f1")      # array of ensemble members
res="1"                   # resolution (degrees) of interpolated grid
remapMethod="bil"         # interpolation method (cdo remap*)
overwriteWeights=false    # if true, re-calculate interp. weights
overwriteRemapped=false   # if true, re-interpolate data even if it exists
keepRemappedSiconc=false  # if true, do not delete interpolated sea ice concentration data

# : after flag indicates that the option requires an argument
while getopts m:x:e:r:c:wok flag
do
    case "${flag}" in
        m) model=${OPTARG};;
        x) experiment=${OPTARG};;
        e) members=(${OPTARG//,/ });;
        r) res=${OPTARG};;
        c) remapMethod=${OPTARG};;
        w) overwriteWeights=true;;
        o) overwriteRemapped=true;;
        k) keepRemappedSiconc=true;;
    esac
done

# Directory containing CMIP6 data. Assumes the same directory
# structure as other processing scripts (i.e., that sea ice
# concentration data resides at
# ${baseDir}/${model}/${experiment}/siconc/*.nc):
baseDir="/storage/research/cpom/gb919150/CMIP6"

# Subdirectory names to contain remapping weights and
# interpolated sea ice concentration data, respectively:
weightsDirName="_remapWeights"
swapDirName="_swap"

weightsDir="${baseDir}/${weightsDirName}"
swapDir="/${baseDir}/${swapDirName}"

# Location of python script:
pyScript="${HOME}/phd/process_cmip6_data/sea_ice"
pyScript="${pyScript}/calc_iel.py"

# Common options to pass to python (before script name):
pyOpts="-u -W ignore"

# Print options to console:
echo "------------------------------------------------"
echo "Ice edge latitude: remapping and calculation"
echo "------------------------------------------------"
echo "Model                     : ${model}"
echo "Experiment                : ${experiment}"
echo ""
echo "Remap method              : ${remapMethod}"
echo "Remap resolution          : ${res} deg"
echo "Overwriting weights       : ${overwriteWeights}"
echo "Overwriting remapped files: ${overwriteRemapped}"
echo "Keeping remapped files    : ${keepRemappedSiconc}"
echo ""
echo "./ = ${baseDir}"
echo "Weights files             : ./${weightsDirName}/"
echo "Intermediate files        : ./${swapDirName}/"
echo "------------------------------------------------"

mkdir -p ${swapDir}
mkdir -p ${weightsDir}

siconcDir="${baseDir}/${model}/${experiment}/siconc"

# Generate weights from first ensemble member, first file:
siconcFiles0=( $(ls ${siconcDir}/siconc_SImon_${model}_${experiment}_${members[0]}_*.nc) )

if [ ! -f "${weightsDir}/${model}_${remapMethod}_weights_${res}deg.nc" ] || [ "${overwriteWeights}" = true ]
then
    echo "Generating weights..."
    cdo gen${remapMethod},global_${res} ${siconcFiles0[0]} \
        ${weightsDir}/${model}_${remapMethod}_weights_${res}deg.nc
else
    echo "Weights already generated:"
    echo "./${weightsDirName}/${model}_${remapMethod}_weights_${res}deg.nc"
fi


for m in ${members[@]}
do
    # Raw files for this ensemble member:
    siconcFilesM=( $(ls ${siconcDir}/siconc_SImon_${model}_${experiment}_${m}_*.nc) )
    
    # Remap each file:
    for x in ${siconcFilesM[@]}
    do
        if [ ! -f "${swapDir}/remapped_${remapMethod}_${res}deg_${x##*/}" ] || [ "${overwriteRemapped}" = true ]
        then
            cdo remap,global_${res},${weightsDir}/${model}_${remapMethod}_weights_${res}deg.nc \
                "${x}" "${swapDir}/remapped_${remapMethod}_${res}deg_${x##*/}"
        else
            echo "Already remapped: ${x##*/}"
        fi
    done
done


python ${pyOpts} ${pyScript} -x ${experiment} -m ${model} \
    --datadir ${swapDir} -g ${res} -r ${remapMethod}


if [ ! "${keepRemappedSiconc}" = true ]
then
    echo "Removing ./${swapDirName}/remapped_${remapMethod}_${res}deg_siconc*${model}_${experiment}_*.nc"
    rm ${swapDir}/remapped_${remapMethod}_${res}deg_siconc*${model}_${experiment}_*.nc
fi


exit 0
