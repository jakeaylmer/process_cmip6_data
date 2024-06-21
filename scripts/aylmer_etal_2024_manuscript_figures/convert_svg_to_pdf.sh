#!/bin/bash

# Convert SVG output from matplotlib to PDF using Inkscape.
# 
# Fig. S4 requires manual editing in the Inkscape GUI
# (specifically, rotation through 180 degrees) before exporting
# to PDF.
# 
# Why not just save directly to PDF from matplotlib in the
# first place? Because there are bugs or limitations with
# matplotlib saving to PDF, specifically with regards to fonts
# (some special characters are converted to paths rather than
# embedded as text*) and the figure layout (it does not seem to
# precisely match the interactively-displayed figure which is
# used to test the layout). Exporting initially to SVG, and
# then using Inkscape to save PDF copies, avoids these issues.
# 
# *known bug: github.com/matplotlib/matplotlib/issues/21797
# ----------------------------------------------------------- # 

figDirIn=./figures
figDirOut=./figures

for fig in "1" "2" "3" "4" "5" "6" "7" "8" "S1" "S2" "S3" "S5" "S6"
do
    if [[ -f ${figDirIn}/fig${fig}.svg ]]
    then
        inkscape ${figDirIn}/fig${fig}.svg \
            --export-type=pdf \
            --export-filename=${figDirOut}/fig${fig}.pdf
	echo "Saved ${figDirIn}/fig${fig}.pdf"
    else
        echo "Figure ${figDirIn}/fig${fig}.svg does not exist; skipping."
    fi
done

exit 0
