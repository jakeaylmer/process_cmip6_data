#!/bin/bash

# Default command-line arguments:
model="AWI-CM-1-1-MR"
experiment="piControl"
members=("r1i1p1f1")
res="0.25"
remapMethod="dis"
overwriteWeights=false
overwriteRemapped=false
keepRemappedSiconc=false

# Working directory; this MUST be set to the directory of the
# package but can be passed as an option (-d):
workDir="${HOME}"

# : after flag indicates that the option requires an argument
while getopts x:e:r:c:d:wok flag
do
    case "${flag}" in
        x) experiment=${OPTARG};;
        e) members=(${OPTARG//,/ });;
        r) res=${OPTARG};;
        c) remapMethod=${OPTARG};;
        d) remapMethod=${OPTARG};;
        w) overwriteWeights=true;;
        o) overwriteRemapped=true;;
        k) keepRemappedSiconc=true;;
    esac
done

cd ${workDir}

# Directory containing CMIP6 data. Assumes the same directory
# structure as other processing scripts (i.e., that sea ice
# concentration data resides at
# ${baseDir}/${model}/${experiment}/siconc/*.nc):
baseDir=$(head -n 1 ./paths/path_cmip6_raw_data.txt)

# Location of python script:
pyScript="${workDir}/sea_ice/calc_iel.py"

weightsDirName="_remapWeights"
swapDirName="_swap"

weightsDir="${baseDir}/${weightsDirName}"
swapDir="/${baseDir}/${swapDirName}"

pyOpts="-u -W ignore"

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


echo "Generating land mask for remapped (${res} deg ${remapMethod}) data"
cdo -f nc -remap${remapMethod},global_${res} -topo ${swapDir}/topo_${res}deg.nc
#cdo -expr,'topo = ((topo<0.0)) ? 1.0 : 0.0' ${swapDir}/topo_${res}deg.nc ${swapDir}/lm_${res}deg.nc
cdo -ltc,0 ${swapDir}/topo_${res}deg.nc ${swapDir}/lm_${res}deg.nc


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
        
        cdo -div "${swapDir}/remapped_${remapMethod}_${res}deg_${x##*/}" \
            ${swapDir}/lm_${res}deg.nc \
            "${swapDir}/masked_remapped_${remapMethod}_${res}deg_${x##*/}"
        
        rm "${swapDir}/remapped_${remapMethod}_${res}deg_${x##*/}"
        mv "${swapDir}/masked_remapped_${remapMethod}_${res}deg_${x##*/}" \
            "${swapDir}/remapped_${remapMethod}_${res}deg_${x##*/}"
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
