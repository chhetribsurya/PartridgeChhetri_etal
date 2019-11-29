#!/bin/bash 

export RUN_PATH=`pwd`
export REFERENCE_DIR="/gpfs/gpfs1/home/schhetri/for_encode/spp/reference"
export GENOME="/gpfs/gpfs1/home/schhetri/for_encode/hg19_genome/hg19-male.fa"
export INPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
export OUTPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/unique_TFs"

### Set BSUB parameters:
export CORE_NUM=2
export MEME_CORES=2

#export BSUB_OPTIONS="-We 24:00 -q c7normal -n $CORE_NUM -R span[hosts=1] -R rusage[mem=25000]" # Using new cluster
export BSUB_OPTIONS="-We 24:00 -q c7normal -R rusage[mem=10000]" # Using new cluster

### Set the path for all the tools to be consistent:
export MEME_SUITE_PATH="/gpfs/gpfs2/software/meme-4.11.3/bin"
export BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin"
export MOTIF_DB_PATH="/gpfs/gpfs2/software/meme-4.11.3/motif_databases"

### generate the null sequences:
export NULL_GENERATE_SCRIPT="$RUN_PATH/nullseq_generate.py"
export NULL_PARAMETERS="-x 2 -r 1 "
export NULL_HG19_INDICES="$RUN_PATH/nullseq_indices_hg19/"

if [[ ! -d $OUTPUT_DIR ]]; then 
	mkdir -p $OUTPUT_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create base output dir: $OUTPUT_DIR. Exiting"; exit 1; fi
fi

export LOG_DIR="$OUTPUT_DIR/log_files"

if [[ ! -d $LOG_DIR ]]; then 
	mkdir -p $LOG_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create log dir: $LOG_DIR. Exiting"; exit 1; fi
fi

for each_idr_file in $(ls $INPUT_DIR/SL*); do
    export TF_NAME=$(basename $each_idr_file | awk 'BEGIN{FS=".";"_";OFS="_"} {print $1,$5,$NF}')
    echo "processing  $TF_NAME ..."

    if [[ ! -d $OUTPUT_DIR/$TF_NAME ]]; then
        mkdir -p $OUTPUT_DIR/$TF_NAME
    fi    

    JOB_NAME=$(echo $TF_NAME | cut -f1,4 -d "_") 
    bsub $BSUB_OPTIONS -J "Comprehensive motif analysis for IDR_${JOB_NAME}" -o $LOG_DIR/comprehensive_motif_calling.out $RUN_PATH/motif_analysis_final_pipeline.sh $each_idr_file $OUTPUT_DIR $TF_NAME

done



