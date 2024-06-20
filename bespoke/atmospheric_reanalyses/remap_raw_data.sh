#!/bin/bash

# This script is used for ERA5 and MERRA-2 which provide
# temperatures at points including exactly 90 degrees_north. So
# that we can use grid cell areas for area integrals in a
# consistent way with the other reanalyses (CFSR, CFSv2, and
# JRA-55), linearly interpolate these data to regular grids of
# the same resolution as the raw data. This makes very little
# difference to the resulting area integrals but puts points at
# 'grid cell centers'.

#rean="ERA5"
#newres="r1440x720"
rean="MERRA-2"
newRes="r576x360"

y1=1980
y2=2023

# Paths to raw and output data:
rawDataDir="/storage/research/cpom/gb919150"
rawDataDir="${rawDataDir}/monthly_near_surface_air_temperature"
rawDataDir="${rawDataDir}/${rean}_raw"

outDataDir="/storage/research/cpom/gb919150"
outDataDir="${outDataDir}/monthly_near_surface_air_temperature"
outDataDir="${outDataDir}/${rean}"

mkdir -p ${outDataDir}

for y in $(seq $y1 $y2)
do
    cdo remapbil,${newRes} \
        ${rawDataDir}/t2.global.monthly.${y}.nc \
        ${outDataDir}/t2.global.monthly.${y}.nc
done

exit 0
