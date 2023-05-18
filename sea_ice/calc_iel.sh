#!/bin/bash

# Default command-line arguments:
model="CESM2-FV2"
experiment="piControl"
res="2"
remapMethod="bil"
overwriteWeights=false
overwriteRemapped=false
keepRemappedSiconc=false

# : after flag indicates that the option requires an argument
while getopts m:x:r:c:wok flag
do
    case "${flag}" in
        m) model=${OPTARG};;
        x) experiment=${OPTARG};;
        r) res=${OPTARG};;
        c) remapMethod=${OPTARG};;
        w) overwriteWeights=true;;
        o) overwriteRemapped=true;;
        k) keepRemappedSiconc=true;;
    esac
done

members=("r1i1p1f1")

baseDir="/storage/basic/cpom/gb919150/CMIP6"

weightsDirName="_remapWeights"
swapDirName="_swap"

weightsDir="${baseDir}/${weightsDirName}"
swapDir="/${baseDir}/${swapDirName}"

pyScript="/home/users/gb919150/phd/process_cmip6_data/sea_ice"
pyScript="${pyScript}/calc_iel.py"
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
