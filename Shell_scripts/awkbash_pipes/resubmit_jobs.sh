#!/bin/bash 

export RUN_PATH=`pwd`
export REFERENCE_DIR="/gpfs/gpfs1/home/schhetri/for_encode/spp/reference"
export GENOME="/gpfs/gpfs1/home/schhetri/for_encode/hg19_genome/hg19-male.fa"
#export INPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/dups"
export INPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
#export OUTPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis/dups"
export OUTPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total_allpeaks/unique_TFs"

### Set BSUB parameters:
#export CORE_NUM=2
#export MEME_CORES=2

#export BSUB_OPTIONS="-We 24:00 -q priority -n $CORE_NUM -R span[hosts=1] -R rusage[mem=25000]" # Using new cluster
#export BSUB_OPTIONS="-We 24:00 -q c7normal -n $CORE_NUM -R span[hosts=1] -R rusage[mem=25000]" # Using new cluster
export BSUB_OPTIONS="-We 24:00 -q c7normal -R rusage[mem=12288]" # Using new cluster
#export JOB_PREFIX="test_batch_I"

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

while read job_info; do
    each_idr_file=$(echo $job_info| awk '{print $2}')
    export TF_NAME=$(basename $each_idr_file | awk 'BEGIN{FS=".";"_";OFS="_"} {print $1,$5,$NF}')
    echo "processing  $TF_NAME ..."

    if [[ ! -d $OUTPUT_DIR/$TF_NAME ]]; then
        mkdir -p $OUTPUT_DIR/$TF_NAME
    fi    

    JOB_NAME=$(echo $TF_NAME | cut -f1,4 -d "_") 
    TF=$(basename $each_idr_file | awk 'BEGIN{FS=".";"_";OFS="_"} {print $NF}' | sed -e 's/\[\([a-zA-Z0-9]*\)\]/_\1/g' )
    #echo $TF
    
    bsub $BSUB_OPTIONS -J "Comp peaks motif analysis for $TF   :::   IDR_${JOB_NAME}" -o $LOG_DIR/resubmit_meme_chip_motif_calling_uniq_TFs.out $RUN_PATH/meme_chip_peaks_motif_analysis_final_allpeaks.sh $each_idr_file $OUTPUT_DIR $TF_NAME
    #bsub -We 24:00 -q priority -n 1 -R span[hosts=1] -J "Call motif analysis with meme_chip" -o ./${JOB_PREFIX}_motif_calling_all_peak.out $RUN_PATH/submit_meme_chip_peaks_motif_analysis_final.sh

done < ./resubmit.txt



