#!/bin/bash

##############################################################
### Run Motif Analysis                                     ###
### Uses Meme scripts and data base                        ###
### Dependent on SPP peak calling finishing sucessfully   ###
### Output is a directory with eps images and html display ###
##############################################################

#SOURCE_DIR="/opt/HAIB/myerslab/etc"
SOURCE_DIR="/gpfs/gpfs2/software/HAIB/myerslab/etc"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then 
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

### Verify we are not running on the head node. ###
if [ -z "$LSB_JOBID" ]; then log_msg fatal "Please run on a compute node. Exiting"; exit 1; fi

#Set the variables:
IDR_PEAK_FILE=$1
OUTPUT_DIR=$2
TF_NAME=$3

#R1=$2
#R2=$3
#LIBRARY=IDR_${R1}_${R2}

### Check if the IDR passed peak file exists, and is more than 100 peaks; but ideally, at least 500 peaks would be needed to call the motifs from those sequences

if [[ -f $IDR_PEAK_FILE ]];then
	wc_num=$(cat $IDR_PEAK_FILE | wc -l)
	peak_num=$(echo $wc_num | cut -d " " -f 1)
	if [[ $peak_num -lt 100 ]]; then
		echo -e "\nToo few IDR passed peaks, so the model cannot be built, skipping motif finding\n"
		PEAK_CALLING_FAILED="TRUE"
	fi
fi

#TEMP_IDR_DIR=$BASE_DIR/IDR_${R1}_${R2}/meme_chip/temp
#MOTIF_OUTPUT=$BASE_DIR/IDR_${R1}_${R2}/meme_chip
TEMP_IDR_DIR=$OUTPUT_DIR/$TF_NAME/meme_chip/temp
MOTIF_OUTPUT=$OUTPUT_DIR/$TF_NAME/meme_chip

if [[ ! -d $TEMP_IDR_DIR ]]; then mkdir -p $TEMP_IDR_DIR; fi
if [[ ! -d $MOTIF_OUTPUT ]]; then mkdir -p $MOTIF_OUTPUT; fi


if [ -z $PEAK_CALLING_FAILED ]; then

    SORTED_FILE=$TEMP_IDR_DIR/${TF_NAME}.sorted_file
    ORIGINAL_FASTA_FILE=$TEMP_IDR_DIR/${TF_NAME}_sorted_file.fasta
    ORIGINAL_BG_FILE=$TEMP_IDR_DIR/${TF_NAME}_sorted_file_background.bed
    ORIGINAL_BG_FASTA_FILE=$TEMP_IDR_DIR/${TF_NAME}_sorted_background_file.fasta
    ORIGINAL_BG_MARKOV_MODEL_FILE=$TEMP_IDR_DIR/${TF_NAME}_sorted_background_file.model
    ORIGINAL_BG_ZERO_MARKOV_MODEL_FILE=$TEMP_IDR_DIR/${TF_NAME}_sorted_background_file.zero_model

    CENTERED_FILE=$TEMP_IDR_DIR/${TF_NAME}_summit_500bp.bed
    CENTERED_FASTA_FILE=$TEMP_IDR_DIR/${TF_NAME}_summit_500bp_file.fasta
    BG_FILE=$TEMP_IDR_DIR/${TF_NAME}_background.bed
    BG_FASTA_FILE=$TEMP_IDR_DIR/${TF_NAME}_background_file.fasta
    BG_MARKOV_MODEL_FILE=$TEMP_IDR_DIR/${TF_NAME}_background_file.model
    BG_ZERO_MARKOV_MODEL_FILE=$TEMP_IDR_DIR/${TF_NAME}_background_file.zero_model

    CENTRIMO_DIR=$OUTPUT_DIR/$TF_NAME/meme_cetrimo_dir
    TOMTOM_DIR=$OUTPUT_DIR/$TF_NAME/meme_tomtom_dir
    FIMO_DIR=$OUTPUT_DIR/$TF_NAME/meme_fimo_dir
    SPAMO_DIR=$OUTPUT_DIR/$TF_NAME/meme_spamo_dir
    MEME_DB_PATH=$MOTIF_OUTPUT/meme_out
    
    if [[ ! -d $CENTRIMO_DIR ]];then mkdir -p $CENTRIMO_DIR; fi
    if [[ ! -d $TOMTOM_DIR ]];then mkdir -p $TOMTOM_DIR; fi
    if [[ ! -d $FIMO_DIR ]];then mkdir -p $FIMO_DIR; fi
    if [[ ! -d $SPAMO_DIR ]];then mkdir -p $SPAMO_DIR; fi

#    ### sorting of the peaks based on signal float value:
#    cat $IDR_PEAK_FILE | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5,$7,$10}' | sort -k 5,5gr > $SORTED_FILE
#    echo -e "\nsorting of IDR peak file completed...\n"   
#
    ### fetching the DNA sequences of sorted regions using hg19 male fasta reference seq:
    $BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $SORTED_FILE -fo $ORIGINAL_FASTA_FILE
    echo -e "\nFasta extraction of sorted bed regions completed...\n"   

    ### Generate 2X null sequences or random sequences with matched GC content, repeat fraction with user input sequence length for sorted file:
    python $NULL_GENERATE_SCRIPT $NULL_PARAMETERS -o $ORIGINAL_BG_FILE $SORTED_FILE hg19 $NULL_HG19_INDICES
    echo -e  "\nGeneration of randomic genomic regions matching GC% and length for sorted bed regions completed...\n"   

    ### fetching the DNA sequences of randomly generated null sequences with matched GC content and repeat fraction using hg19 ref seq:
    $BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $ORIGINAL_BG_FILE -fo $ORIGINAL_BG_FASTA_FILE
    echo -e "\nFasta extraction of sorted original background file completed...\n"   

    ### create fasta markov model as a background file:
    $MEME_SUITE_PATH/fasta-get-markov -m 1 $ORIGINAL_BG_FASTA_FILE $ORIGINAL_BG_MARKOV_MODEL_FILE
    echo -e "\nMarkov model generation for sorted original background file completed...\n"   

    ### create zero order fasta markov model as a background file (useful for spamo analysis, if motif spacing analysis performed on original file rather than centered ones):
    $MEME_SUITE_PATH/fasta-get-markov -m 0 $ORIGINAL_BG_FASTA_FILE $ORIGINAL_BG_ZERO_MARKOV_MODEL_FILE
    echo -e "\nZero order Markov model generation for sorted original background file completed...\n"   

#    ### generating bed file of 500bp regions centered on peak-summits: 
#    awk 'BEGIN{OFS="\t"} {chromStart=$2; summit=$6; midPos=chromStart+summit; print $1, midPos-250, midPos+250;}' $SORTED_FILE > $CENTERED_FILE
#    echo -e "\nCentering of IDR peak file 250 upstream and 250 downstream completed...\n"   
#    
#    ### Generate 2X null sequences or random sequences with matched GC content, repeat fraction with user input sequence length:
#    python $NULL_GENERATE_SCRIPT $NULL_PARAMETERS -o $BG_FILE $CENTERED_FILE hg19 $NULL_HG19_INDICES
#    echo -e  "\nGeneration of randomic genomic regions matching GC% and length completed...\n"   
#
#    ### fetching the DNA sequences of peak summit regions using hg19 male fasta reference seq:
#    $BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $CENTERED_FILE -fo $CENTERED_FASTA_FILE
#    echo -e "\nFasta extraction of summit centered file completed...\n"   
#
#    ### fetching the DNA sequences of randomly generated null sequences with matched GC content and repeat fraction using hg19 ref seq:
#    $BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $BG_FILE -fo $BG_FASTA_FILE
#    echo -e "\nFasta extraction of background file completed...\n"   
#
#    ### create fasta markov model as a background file:
#    $MEME_SUITE_PATH/fasta-get-markov -m 1 $BG_FASTA_FILE $BG_MARKOV_MODEL_FILE
#    echo -e "\nMarkov model generation for background file completed...\n"   
#
    ### create zero order fasta markov model as a background file, useful esp. for spamo:
    $MEME_SUITE_PATH/fasta-get-markov -m 0 $BG_FASTA_FILE $BG_ZERO_MARKOV_MODEL_FILE
    echo -e "\nZero order Markov model generation for background file completed...\n"   
#
#    ### Run meme-chip with parallel cores and nonrandom features since the input fasta is sorted with decreasing confidence of the peaks or signal float value:
#    echo -e "\nStarting the meme-chip operation using the centered fasta file and bg markov model file generated in prior steps...\n"
#
#    $MEME_SUITE_PATH/meme-chip -oc $MOTIF_OUTPUT -index-name meme_combined.html -db $MOTIF_DB_PATH/JASPAR/JASPAR_CORE_2016_vertebrates.meme -db $MOTIF_DB_PATH/MOUSE/uniprobe_mouse.meme -dna -bfile $BG_MARKOV_MODEL_FILE -norand -nmeme 500 -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 5 -meme-maxsize 100000 -dreme-m 5 -centrimo-local -centrimo-score 5 -centrimo-ethresh 10 ${CENTERED_FASTA_FILE}
#    echo -e "\nMeme-chip analysis completed....\n"

#    $MEME_SUITE_PATH/centrimo --seqlen 500 --verbosity 1 --oc $CENTRIMO_DIR --neg $BG_FASTA_FILE --disc --bfile $BG_MARKOV_MODEL_FILE --local --score 5 --ethresh 10 $CENTERED_FASTA_FILE $MEME_DB_PATH/meme.xml $MOTIF_DB_PATH/JASPAR/JASPAR_CORE_2016_vertebrates.meme $MOTIF_DB_PATH/MOUSE/uniprobe_mouse.meme
#    echo -e "\nCentrimo analysis completed...\n"

    $MEME_SUITE_PATH/tomtom -verbosity 1 -oc $TOMTOM_DIR -png -eps -min-overlap 5 -dist pearson -evalue -thresh 1 -no-ssc -bfile $BG_MARKOV_MODEL_FILE $MEME_DB_PATH/meme.xml $MOTIF_DB_PATH/JASPAR/JASPAR_CORE_2016_vertebrates.meme $MOTIF_DB_PATH/MOUSE/uniprobe_mouse.meme 
    echo -e "\nTOMTOM analysis completed...\n"

    for meme_id in {1..5}; do
        SPAMO_OUTDIR=${SPAMO_DIR}/spamo_motif_${meme_id}
        if [[ ! -d ${SPAMO_OUTDIR} ]]; then mkdir -p ${SPAMO_OUTDIR}; fi
        echo -e "\nProcessing Spamo analysis on MEME $meme_id ...\n"
        ### If interested in scanning the motif with the original file rather than centered region(but centered region should be better inorder to be consistent with the range scanning of sec motif from primary motif):
        #$MEME_SUITE_PATH/spamo -verbosity 1 -oc ${SPAMO_OUTDIR} -range 100 -png -eps -bgfile $ORIGINAL_BG_ZERO_MARKOV_MODEL_FILE -primary $meme_id ${ORIGINAL_FASTA_FILE} $MEME_DB_PATH/meme.xml $MEME_DB_PATH/meme.xml $MOTIF_DB_PATH/JASPAR/JASPAR_CORE_2016_vertebrates.meme $MOTIF_DB_PATH/MOUSE/uniprobe_mouse.meme
        $MEME_SUITE_PATH/spamo -verbosity 1 -oc ${SPAMO_OUTDIR} -range 100 -png -eps -bgfile $BG_ZERO_MARKOV_MODEL_FILE -primary $meme_id ${CENTERED_FASTA_FILE} $MEME_DB_PATH/meme.xml $MEME_DB_PATH/meme.xml $MOTIF_DB_PATH/JASPAR/JASPAR_CORE_2016_vertebrates.meme $MOTIF_DB_PATH/MOUSE/uniprobe_mouse.meme
    done
    echo -e "\nSpamo analysis completed...\n"

    ### Make sure that you scan the motif within the original sequence and not on the 500bp centered fasta file, which is useful for centrimo and spamo purpose only:
    #$MEME_SUITE_PATH/fimo --parse-genomic-coord --verbosity 1 --oc ${FIMO_DIR}_w_ZeroOrder --bgfile $ORIGINAL_BG_ZERO_MARKOV_MODEL_FILE $MEME_DB_PATH/meme.xml ${ORIGINAL_FASTA_FILE} ##for scannning with zero order bg file
    $MEME_SUITE_PATH/fimo --parse-genomic-coord --verbosity 1 --oc $FIMO_DIR --bgfile $ORIGINAL_BG_MARKOV_MODEL_FILE $MEME_DB_PATH/meme.xml ${ORIGINAL_FASTA_FILE} 
    echo -e "\nFIMO analysis completed...\n"

fi

echo -e "\nMotif analysis completed!!!\n"




