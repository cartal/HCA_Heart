#!/bin/bash -x

for i in $(cat samples.txt)

do

    bash ct5_cellranger_pipeline.sh cram <REFERENCE> sample $i

done

