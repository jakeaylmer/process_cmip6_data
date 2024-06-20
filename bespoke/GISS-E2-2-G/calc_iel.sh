#!/bin/bash

# Default command-line arguments:
model="GISS-E2-2-G"
experiment="piControl"
members=("r1i1p3f1")
remapMethod="dis"
overwriteLandMask=false
overwriteMaskedSiconc=false
keepLandMask=false
keepMaskedSiconc=false

# Working directory; this MUST be set to the directory of the
# package but can be passed as an option (-d):
workDir="${HOME}"

# : after flag indicates that the option requires an argument
while getopts x:e:c:d:klow flag
do
    case "${flag}" in
        x) experiment=${OPTARG};;
        e) members=(${OPTARG//,/ });;
        c) remapMethod=${OPTARG};;
        d) workDir=${OPTARG};;
        k) keepMaskedSiconc=true;;
        l) keepLandMask=true;;
        o) overwriteLandMask=true;;
        w) overwriteMaskedSiconc=true;;
    esac
done

# Directory containing CMIP6 data. Assumes the same directory
# structure as other processing scripts (i.e., that sea ice
# concentration data resides at
# ${baseDir}/${model}/${experiment}/siconc/*.nc):
baseDir=$(head -n 1 ./paths/path_cmip6_raw_data.txt)

# Location of python script:
pyScript="${workDir}/bespoke/GISS-E2-2-G/calc_iel.py"

swapDirName="_swap"
swapDir="${baseDir}/${swapDirName}"

pyOpts="-u -W ignore"

echo "------------------------------------------------"
echo "Ice edge latitude: remapping and calculation"
echo "------------------------------------------------"
echo "Model                    : ${model}"
echo "Experiment               : ${experiment}"
echo ""
echo "Remap method (land mask) : ${remapMethod}"
echo "Keeping land-mask file   : ${keepLandMask}"
echo "Keeping masked files     : ${keepMaskedSiconc}"
echo "Overwriting land mask    : ${overwriteLandMask}"
echo "Overwriting masked siconc: ${overwriteMaskedSiconc}"
echo ""
echo "./ = ${baseDir}"
echo "Intermediate files       : ./${swapDirName}/"
echo "------------------------------------------------"

mkdir -p ${swapDir}

siconcDir="${baseDir}/${model}/${experiment}/siconca"
siconcFiles0=( $(ls ${siconcDir}/siconca_SImon_${model}_${experiment}_${members[0]}_*.nc) )

if [ ! -f "${swapDir}/lmask_${remapMethod}.nc" ] || [ "${overwriteLandMask}" = true ]
then
    echo "Generating land mask (using remap${remapMethod})"
    cdo -f nc -remap${remapMethod},${siconcFiles0[0]} -topo ${swapDir}/topo.nc
    cdo -ltc,0 ${swapDir}/topo.nc ${swapDir}/lmask_${remapMethod}.nc
    rm ${swapDir}/topo.nc
else
    echo "Already exists: ${swapDir}/lmask_${remapMethod}.nc"
fi

for m in ${members[@]}
do
    # Raw files for this ensemble member:
    siconcFilesM=( $(ls ${siconcDir}/siconca_SImon_${model}_${experiment}_${m}_*.nc) )
    
    # Remap each file:
    for x in ${siconcFilesM[@]}
    do
        if [ ! -f "${swapDir}/masked_${x##*/}" ] || [ "${overwriteMaskedSiconc}" = true ]
        then
            cdo -div ${siconcDir}/"${x##*/}" \
                ${swapDir}/lmask_${remapMethod}.nc \
                "${swapDir}/masked_${x##*/}"
        else
            echo "Already exists: ${swapDir}/masked_${x##*/}"
        fi
        
    done
done


python ${pyOpts} ${pyScript} -x ${experiment} -m ${model} \
    --datadir ${swapDir}

if [ ! "${keepMaskedSiconc}" = true ]
then
    echo "Removing ./${swapDirName}/masked_siconca*${model}_${experiment}_*.nc"
    rm ${swapDir}/masked_siconca*${model}_${experiment}_*.nc
fi

if [ ! "${keepLandMask}" = true ]
then
    echo "Removing ./${swapDirName}/lmask_${remapMethod}.nc"
    rm ${swapDir}/lmask_${remapMethod}.nc
fi

exit 0
