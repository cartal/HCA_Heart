#!/bin/bash -x

### SET UP ENVIRONMENT

E=200
INPUT_DIR=/mnt/data/hca_heart/cellranger/harvard/premrna/
OUTPUT_DIR=/mnt/data/hca_heart/cellbender/harvard/premrna/
REPORTS=/mnt/data/hca_heart/cellranger/harvard/premrna/reports/

### RUN CELLBENDER
# For this to work, you will need a file with at least three columns:
#sampleID,expected_cells,included_cells


for i in $(echo $(cat samples.txt) | cut -f1 -d-)

do

    EXPECTED=$(echo $(cat samples.txt) | cut -f2 -d-)
    INCLUDED=$(echo $(cat samples.txt) | cut -f3 -d-)

    ln -s $INPUT_DIR/$i/outs/$i.h5
    
    cellbender remove-background --input $i.h5 --output $i.cellbender.h5 --expected-cells $EXPECTED --total-droplets-included $INCLUDED --epochs $E --cuda

    rm $i.h5

    mv *.cellbender.h5 $OUTPUT_DIR

    mv *.log *.pdf $REPORTS

done
