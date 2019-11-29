#!/bin/bash 

export RUN_PATH=`pwd`
export REFERENCE_DIR="/gpfs/gpfs1/home/schhetri/for_encode/spp/reference"
export GENOME="/gpfs/gpfs1/home/schhetri/for_encode/hg19_genome/hg19-male.fa"
export INPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/idr_passed_peaks_total/unique_TFs"
export OUTPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/finalEncode_reduced_cisbp_jasparbased_unique_TFanalysis"

### Set BSUB parameters:
export CORE_NUM=2
export MEME_CORES=2

#export BSUB_OPTIONS="-We 24:00 -q c7normal -n $CORE_NUM -R span[hosts=1] -R rusage[mem=25000]" # Using new cluster
#export BSUB_OPTIONS="-We 24:00 -q c7normal -R rusage[mem=10000]" # Using new cluster
export BSUB_OPTIONS="-We 24:00 -q c7normal" # Using new cluster

### Set the path for all the tools to be consistent:
#export MEME_SUITE_PATH="/gpfs/gpfs2/software/meme-4.11.3/bin"
export MEME_SUITE_PATH="/gpfs/gpfs2/software/meme-4.11.4/bin"
export BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin"
export MOTIF_DB_PATH="/gpfs/gpfs1/home/schhetri/Tools/final_custom_CisBP_database"

### generate the null sequences:
export NULL_GENERATE_SCRIPT="$RUN_PATH/nullseq_generate.py"
export NULL_PARAMETERS="-x 2 -r 1 "
export NULL_HG19_INDICES="$RUN_PATH/nullseq_indices_hg19/"

if [[ ! -d $OUTPUT_DIR ]]; then 
	mkdir -p $OUTPUT_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create base output dir: $OUTPUT_DIR. Exiting"; exit 1; fi
fi

export LOG_DIR="$OUTPUT_DIR/log_files_stop"

if [[ ! -d $LOG_DIR ]]; then 
	mkdir -p $LOG_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create log dir: $LOG_DIR. Exiting"; exit 1; fi
fi

for each_idr_file in $(ls $INPUT_DIR/SL*); do
    export TF_NAME=$(basename $each_idr_file | awk 'BEGIN{FS=".";"_";OFS="_"} {print $1,$5,$NF}')
    export TF_NAME_IDENTIFIER=$(echo $TF_NAME | cut -f 6 -d "_"|sed 's/\[FLAG\]//')
    echo "processing ${TF_NAME_IDENTIFIER} : $TF_NAME ..."

    if [[ ! -d $OUTPUT_DIR/$TF_NAME ]]; then
        mkdir -p $OUTPUT_DIR/$TF_NAME
    fi    

    JOB_NAME=$(echo $TF_NAME | cut -f1,4 -d "_") 
    bsub $BSUB_OPTIONS -J "motif similarity analysis for ${TF_NAME_IDENTIFIER}_IDR" -o $LOG_DIR/comp_cisbp_allmotifs_tomtom_calling_stop_Eval100.out $RUN_PATH/cisbpjasp_motif_analysis_pipe.sh $each_idr_file ${OUTPUT_DIR} ${TF_NAME} ${TF_NAME_IDENTIFIER}

done



