#!/bin/bash 

###IDR program requires python3, so new environment created using conda -create, 
###given anaconda2.7 is already created:
#conda create -n python3 python=3.7 anaconda

#First activate anaconda version of python
#source activate python3

export RUN_PATH=`pwd`
export TF_FILE_PATH="/gpfs/gpfs1/home/schhetri/k562_chip/SE_libraries_info"
#export LIB_PATH="/gpfs/gpfs1/home/schhetri/for_encode/spp/Libraries"
export LIB_PATH="/gpfs/gpfs1/home/schhetri/k562_chip/k562_SE_fastq_libraries"
export BASE_DIR="/gpfs/gpfs1/home/schhetri/k562_chip/k562_SE_chip_analysis"
export REFERENCE_DIR="/gpfs/gpfs1/home/schhetri/for_encode/spp/reference"
export GENOME="/gpfs/gpfs1/home/schhetri/for_encode/hg19_genome/hg19-male.fa"
export TEMP_DIR="$BASE_DIR/temp"
export LOG_DIR="$BASE_DIR/log_files"
export BSUB_OPTIONS="-We 24:00 -q c7normal -R rusage[mem=40000]" # Using new cluster
#export BSUB_OPTIONS="-We 24:00 -q normal -n 6 -R span[hosts=1]" # Using new cluster
export JOB_PREFIX="K562_SE"

### Set BSUB parameters:
#export CORE_NUM=2
#export MEME_CORES=2
#export BSUB_OPTIONS="-We 10:00 -n 8 -R span[hosts=1]"
#export BSUB_OPTIONS="-We 10:00 -q night -n 8 -R span[hosts=1]"
#export BSUB_OPTIONS="-We 24:00 -q c7normal -n $CORE_NUM -R span[hosts=1] -R rusage[mem=25000]" # Using new cluster


### Set the path for all the tools to be consistent:
#export MEME_SUITE_PATH="/gpfs/gpfs2/software/meme-4.11.3/bin"
export MEME_SUITE_PATH="/gpfs/gpfs2/software/meme-4.11.4/bin"
export BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin"
export MOTIF_DB_PATH="/gpfs/gpfs2/software/meme-4.11.3/motif_databases"

### generate the null sequences:
export NULL_GENERATE_SCRIPT="$RUN_PATH/nullseq_generate.py"
export NULL_PARAMETERS="-x 2 -r 1 "
export NULL_HG19_INDICES="$RUN_PATH/nullseq_indices_hg19/"

if [[ ! -d $BASE_DIR ]]; then 
	mkdir -p $BASE_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create base output dir: $BASE_DIR. Exiting"; exit 1; fi
fi

if [[ ! -d $TEMP_DIR ]]; then 
	mkdir -p $TEMP_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create temp dir: $TEMP_DIR. Exiting"; exit 1; fi
fi

if [[ ! -d $LOG_DIR ]]; then 
	mkdir -p $LOG_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create log dir: $LOG_DIR. Exiting"; exit 1; fi
fi

if [[ ! -d $REFERENCE_DIR ]]; then 
	mkdir -p $REFERENCE_DIR
	if [[ $? -ne "0" ]]; then echo "Could not create reference dir: $REFERENCE_DIR. Exiting"; exit 1; fi
fi

#export TF_FILE_PATH="/gpfs/gpfs1/home/schhetri/for_encode/spp"
export REP1_FILE_NAME="$TF_FILE_PATH/${JOB_PREFIX}_tf_rep1.txt"
export REP2_FILE_NAME="$TF_FILE_PATH/${JOB_PREFIX}_tf_rep2.txt"
export CONTROL1_FILE_NAME="$TF_FILE_PATH/${JOB_PREFIX}_tf_control1.txt"
export CONTROL2_FILE_NAME="$TF_FILE_PATH/${JOB_PREFIX}_tf_control2.txt"

cat $REP1_FILE_NAME $REP2_FILE_NAME $CONTROL1_FILE_NAME $CONTROL2_FILE_NAME > $TF_FILE_PATH/${JOB_PREFIX}_all_tf.txt
export ALL_TF_FILE_NAME="$TF_FILE_PATH/${JOB_PREFIX}_all_tf.txt"

###########################################

# Only cut and past the below codes within the double hashed linw on the relevant script for reading the array from the file:
#So, to overcome the issue of exporting the bash array variable, so directly: (Read lines in file into an bash array)
readarray LIB_REP1 < $REP1_FILE_NAME
readarray LIB_REP2 < $REP2_FILE_NAME
readarray LIB_CONTROL1 < $CONTROL1_FILE_NAME
readarray LIB_CONTROL2 < $CONTROL2_FILE_NAME
readarray ALL_LIB < $ALL_TF_FILE_NAME

#Reading or geting unique values from an bash array, and then note it's reconverted into bash array by placing array bracket ( ):
LIB_REP1=( $(echo ${LIB_REP1[@]} | tr " " "\n" | tr "\n" " "))
LIB_REP2=( $(echo ${LIB_REP2[@]} | tr " " "\n" | tr "\n" " "))
LIB_CONTROL1=( $(echo ${LIB_CONTROL1[@]} | tr " " "\n" | tr "\n" " "))
LIB_CONTROL2=( $(echo ${LIB_CONTROL2[@]} | tr " " "\n" | tr "\n" " "))
UNIQ_ALL_LIB=( $(echo ${ALL_LIB[@]} | tr " " "\n" | sort -u | tr "\n" " ") )

###########################################

###Sanity check for the different combination of the SL libraries:
for i in "${!LIB_REP1[@]}"; do

	CURRENT_REP1=${LIB_REP1[$i]}
	CURRENT_REP2=${LIB_REP2[$i]}
	CURRENT_CONTROL1=${LIB_CONTROL1[$i]}
	CURRENT_CONTROL2=${LIB_CONTROL2[$i]}
	echo -e "Processing ${CURRENT_REP1}_${CURRENT_REP2}_${CURRENT_CONTROL1}_${CURRENT_CONTROL2} library combination\n"
	
done

bsub -We 24:00 -n 2 -R span[hosts=1] -J "step1_tag_align" -o $LOG_DIR/${JOB_PREFIX}_tag_align.out $RUN_PATH/submit_fastq2tagalign.sh

bsub -We 24:00 -q c7normal -n 1 -R span[hosts=1] -J "step2_pool_peak_call" -o $LOG_DIR/${JOB_PREFIX}_pooling_and_peak_calling.out $RUN_PATH/submit_peak_calls.sh

#bsub -We 24:00 -n 1 -R span[hosts=1] -J "Call IDR calling" -o $LOG_DIR/${JOB_PREFIX}_idr_calling.out $RUN_PATH/submit_run_idr.sh

bsub -We 24:00 -n 1 -R span[hosts=1] -J "Call IDR calling 0.5" -o $LOG_DIR/${JOB_PREFIX}_idr_calling.out $RUN_PATH/submit_run_idr_05.sh

#bsub -We 24:00 -n 1 -R span[hosts=1] -J "Call motif analysis" -o $LOG_DIR/${JOB_PREFIX}_motif_calling.out $RUN_PATH/submit_meme_peaks_motif_analysis.sh

bsub -We 24:00 -q c7normal -n 1 -R span[hosts=1] -J "Call motif analysis with meme_chip" -o ./${JOB_PREFIX}_motif_calling_all_peak.out $RUN_PATH/submit_meme_chip_peaks_motif_analysis_final.sh


