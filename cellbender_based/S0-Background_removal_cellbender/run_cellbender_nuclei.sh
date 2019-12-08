#!/bin/bash -x

### SET UP ENVIRONMENT

E=200
INPUT_DIR=/mnt/data/hca_heart/cellranger/harvard/premrna/
OUTPUT_DIR=/mnt/data/hca_heart/cellbender/harvard/premrna/
REPORTS=/mnt/data/hca_heart/cellranger/harvard/premrna/reports/

### RUN CELLBENDER

# After talking to the developers, their recommendation was to use the number of targeted cells, as 10X recommends 
# and a value within the background plateau that makes sense for the whole expriment (donor). 
# In the case of cells, their recommendation was EXPECTED = 5000, INCLUDED = 13000. 


for i in $(cat samples.txt)

do

    EXPECTED=
    INCLUDED=

    ln -s $INPUT_DIR/$i/outs/$i.h5
    
    cellbender remove-background --input $i.h5 --output $i.cellbender.h5 --expected-cells $EXPECTED --total-droplets-included $INCLUDED --epochs $E --cuda

    rm $i.h5

    mv *.cellbender.h5 $OUTPUT_DIR

    mv *.log *.pdf $REPORTS

done
