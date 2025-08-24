#!/bin/bash

WD=$1
OUTDIR=$2

export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/  /scratch1/fs1/martyomov:/scratch1/fs1/martyomov /home/carisa:/home/carisa" 

mkdir -p $OUTDIR
SAMPLES=(E1_CD8 E1_TOTAL M1_TOTAL M1_CD8 L1_TOTAL L1_CD8)
# 

for SAMPLE in ${SAMPLES[@]}; do
    mkdir -p $OUTDIR/$SAMPLE
    LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q martyomov -G compute-martyomov \
        -J ${SAMPLE}_align -n 8 -M 64GB -o ${SAMPLE}_align.out \
	-e ${SAMPLE}_align.err -R 'select[mem>64MB] rusage[mem=64GB] span[hosts=1]' \
        -a "docker(kalisaz/cellranger:v8.0.1)" /bin/bash -c "cellranger multi --id $SAMPLE \
        --csv=${WD}/alignment/${SAMPLE}_meta.csv \
        --output-dir $OUTDIR/${SAMPLE} \
        --localmem=64 \
        --localcores=16"
done