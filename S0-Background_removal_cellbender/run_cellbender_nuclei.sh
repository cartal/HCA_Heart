#!/bin/bash -x

### SET UP ENVIRONMENT

E=200
INPUT_DIR=/mnt/data/hca_heart/cellranger/harvard/premrna/
OUTPUT_DIR=/mnt/data/hca_heart/cellbender/harvard/premrna/
REPORTS=/mnt/data/hca_heart/cellranger/harvard/premrna/reports/

### RUN CELLBENDER

for i in $(cat samples.txt)

do

    ln -s $INPUT_DIR/$i/outs/$i.h5
    
    cellbender remove-background --input $i.h5 --output $i.cellbender.h5 --expected-cells <INSERT> --total-droplets-included <INSERT> --epochs $E --cuda

    rm $i.h5

    mv *.cellbender.h5 $OUTPUT_DIR

    mv *.log *.pdf $REPORTS

done
