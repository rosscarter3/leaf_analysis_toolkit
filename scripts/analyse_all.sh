#!/usr/bin/bash



for d in ${2}/*; do
    if [[ ${d} == *'BL'* ]];
    then
        echo ${d}

        if [ ${1} == 'segment' ]; then
            python ./ws_segment_from_dir.py ${d}
        elif [ ${1} == 'extract' ]; then
            python ./extract_celldata.py ${d}
        elif [ ${1} == 'heatmap' ]; then
            python ./heatmap.py ${d}
        else
            echo "No process selected"
        fi

#        for f in ${d}/*; do
#            echo ${f}
#            if [[ ${f} == *'seg.png' ]];
#            then
#                rm -r ${f}
#            fi
#        done


    fi
done