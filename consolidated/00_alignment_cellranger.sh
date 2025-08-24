#!/bin/sh

#All paths are assumed to be relative to the workding directory

WD=$1
INDIR=$2
OUTDIR="align_outs/"

GENOME="mm"
export LSF_DOCKER_VOLUMES="/storage1/fs1/martyomov/Active/:/storage1/fs1/martyomov/Active/  /scratch1/fs1/martyomov:/scratch1/fs1/martyomov /home/carisa:/home/carisa" 

WHITELIST="/storage1/fs1/martyomov/Active/IndividualBackUps/carisa/CellRanger_barcodes/barcodes/3M-february-2018.txt"

case "$GENOME" in
	"hs")
		REF_DIR='/storage1/fs1/martyomov/Active/References/10X/SC/Human/refdata-gex-GRCh38-2024-A/'
	;;
	"mm")
		REF_DIR='/storage1/fs1/martyomov/Active/References/10X/SC/Mouse/refdata-gex-GRCm39-2024-A/'
	;;
esac


cd $WD
mkdir -p $OUTDIR



#SAMPLES=( $(ls $INDIR) )
SAMPLES=("KO_HC" "WT_HC")
for SAMPLE in ${SAMPLES[@]}; do
    FASTQ_PATH=${INDIR}/${SAMPLE}/
    mkdir -p $OUTDIR/$SAMPLE/
    LSF_DOCKER_PRESERVE_ENVIRONMENT=false bsub -q martyomov -G compute-martyomov \
        -J ${SAMPLE}_align -n 8 -M 64GB -o ${SAMPLE}_align.out \
	    -e ${SAMPLE}_align.err -R 'select[mem>64MB] rusage[mem=64GB] span[hosts=1]' \
        -a "docker(kalisaz/cellranger:v8.0.1)" /bin/bash -c "cellranger count --id $SAMPLE \
        --output-dir $OUTDIR/$SAMPLE/ \
        --sample $SAMPLE \
        --transcriptome $REF_DIR \
        --fastqs $FASTQ_PATH \
        --create-bam true \
        --localmem=64 \
        --localcores=16" 
        

done