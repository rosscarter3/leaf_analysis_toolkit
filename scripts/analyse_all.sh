#!/usr/bin/bash



for d in ${1}/*; do
    if [[ ${d} == *'BL'* ]];
    then
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

#        for f in ${d}/*; do
#            echo ${f}
#            if [[ ${f} == *'seg.png' ]];
#            then
#                rm -r ${f}
#            fi
#        done


    fi
done