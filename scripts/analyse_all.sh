#!/bin/bash

for d in ${1}/*; do
    if [[ ${d} == *'proj_'* ]]; then
        echo ${d}
        if [ ${2} == 'project' ]; then
            python ./stack2proj.py ${d}
        elif [ ${2} == 'segment' ]; then
            python ./ws_segment_from_dir.py ${d}
        elif [ ${2} == 'extract' ]; then
            python ./extract_celldata.py ${d}
        elif [ ${2} == 'heatmap' ]; then
            python ./heatmap.py ${d}
        else
            echo "No process selected"
        fi
    fi
done
