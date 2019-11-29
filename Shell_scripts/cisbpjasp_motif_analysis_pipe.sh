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
TF_NAME_IDENTIFIER=$4

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

TEMP_IDR_DIR=$OUTPUT_DIR/$TF_NAME/meme_chip/temp
PRIOR_MOTIF_OUTPUT=/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/unique_TFs/$TF_NAME/meme_chip

if [[ ! -d $TEMP_IDR_DIR ]]; then mkdir -p $TEMP_IDR_DIR; fi
#if [[ ! -d $MOTIF_OUTPUT ]]; then mkdir -p $MOTIF_OUTPUT; fi


if [ -z $PEAK_CALLING_FAILED ]; then
    
    # Make sure this file is present there or sort the file, whereas rest of file is generated on the process:
    #SORTED_FILE=$TEMP_IDR_DIR/${TF_NAME}.sorted_file
    #sort -k1,1 -k2,2n ${IDR_PEAK_FILE} > ${SORTED_FILE}

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

    CENTRIMO_DIR=$OUTPUT_DIR/$TF_NAME/meme_centrimo_dir
    CENTRIMO_CISBP_DIR=$OUTPUT_DIR/$TF_NAME/meme_cisbp_centrimo_dir
    CENTRIMO_JASPAR_DIR=$OUTPUT_DIR/$TF_NAME/meme_jaspar_centrimo_dir
    CENTRIMO_CISBP_JASPAR_DIR=$OUTPUT_DIR/$TF_NAME/meme_cisbp_jaspar_centrimo_dir

    TOMTOM_DIR=$OUTPUT_DIR/$TF_NAME/meme_tomtom_dir
    TOMTOM_CISBP_DIR=$OUTPUT_DIR/$TF_NAME/meme_cisbp_tomtom_dir
    TOMTOM_COMP_CISBP_DIR=$OUTPUT_DIR/$TF_NAME/meme_finalredo_cisbp_tomtom_dir_stop
    TOMTOM_JASPAR_DIR=$OUTPUT_DIR/$TF_NAME/meme_jaspar_tomtom_dir_stop
    TOMTOM_JASPAR2018_DIR=$OUTPUT_DIR/$TF_NAME/meme_jaspar_tomtom_dir_2018_stop
    TOMTOM_JASPAR2016_DIR=$OUTPUT_DIR/$TF_NAME/meme_jaspar_tomtom_dir_2016_stop

    FIMO_DIR=$OUTPUT_DIR/$TF_NAME/meme_fimo_dir
    FIMO_CISBP_DIR=$OUTPUT_DIR/$TF_NAME/fimo_cisbp_dir
    FIMO_JASPAR_DIR=$OUTPUT_DIR/$TF_NAME/fimo_jaspar_dir
    #FIMO_CISBP_JASPAR_DIR=$OUTPUT_DIR/$TF_NAME/meme_fimo_cisbp_jaspar_dir

    CISBP_MEME_DIR=$OUTPUT_DIR/$TF_NAME/cisbp_singlefactor_based_meme
    JASPAR_MEME_DIR=$OUTPUT_DIR/$TF_NAME/jaspar_singlefactor_based_meme

    #SPAMO_DIR=$OUTPUT_DIR/$TF_NAME/meme_spamo_dir
    MEME_DBFILE_PATH=$PRIOR_MOTIF_OUTPUT/meme_out
    
    if [[ ! -d $CENTRIMO_DIR ]];then mkdir -p $CENTRIMO_DIR; fi
    if [[ ! -d $CENTRIMO_CISBP_DIR ]];then mkdir -p $CENTRIMO_CISBP_DIR; fi
    if [[ ! -d $CENTRIMO_JASPAR_DIR ]];then mkdir -p $CENTRIMO_JASPAR_DIR; fi
    if [[ ! -d $CENTRIMO_CISBP_JASPAR_DIR ]];then mkdir -p $CENTRIMO_CISBP_JASPAR_DIR; fi

    if [[ ! -d $TOMTOM_DIR ]];then mkdir -p $TOMTOM_DIR; fi
    if [[ ! -d $TOMTOM_CISBP_DIR ]];then mkdir -p $TOMTOM_CISBP_DIR; fi
    if [[ ! -d $TOMTOM_COMP_CISBP_DIR ]];then mkdir -p $TOMTOM_COMP_CISBP_DIR; fi
    if [[ ! -d $TOMTOM_JASPAR_DIR ]];then mkdir -p $TOMTOM_JASPAR_DIR; fi

    if [[ ! -d $FIMO_DIR ]];then mkdir -p $FIMO_DIR; fi
    if [[ ! -d $FIMO_CISBP_DIR ]];then mkdir -p $FIMO_CISBP_DIR; fi
    if [[ ! -d $FIMO_JASPAR_DIR ]];then mkdir -p $FIMO_JASPAR_DIR; fi

    if [[ ! -d $CISBP_MEME_DIR ]];then mkdir -p $CISBP_MEME_DIR; fi
    if [[ ! -d $JASPAR_MEME_DIR ]];then mkdir -p $JASPAR_MEME_DIR; fi
    #if [[ ! -d $SPAMO_DIR ]];then mkdir -p $SPAMO_DIR; fi

##    ### sorting of the peaks based on signal float value:
#    cat $IDR_PEAK_FILE | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$5,$7,$10}' | sort -k 5,5gr > $SORTED_FILE
#    echo -e "\nsorting of IDR peak file completed...\n"   
#
#    ### fetching the DNA sequences of sorted regions using hg19 male fasta reference seq:
#    $BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $SORTED_FILE -fo $ORIGINAL_FASTA_FILE
#    echo -e "\nFasta extraction of sorted bed regions completed...\n"   
#
#    ### Generate 2X null sequences or random sequences with matched GC content, repeat fraction with user input sequence length for sorted file:
#    python $NULL_GENERATE_SCRIPT $NULL_PARAMETERS -o $ORIGINAL_BG_FILE $SORTED_FILE hg19 $NULL_HG19_INDICES
#    echo -e  "\nGeneration of randomic genomic regions matching GC% and length for sorted bed regions completed...\n"   
#
#    ### fetching the DNA sequences of randomly generated null sequences with matched GC content and repeat fraction using hg19 ref seq:
#    $BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $ORIGINAL_BG_FILE -fo $ORIGINAL_BG_FASTA_FILE
#    echo -e "\nFasta extraction of sorted original background file completed...\n"   
#
#    ### create fasta markov model as a background file:
#    $MEME_SUITE_PATH/fasta-get-markov -m 1 $ORIGINAL_BG_FASTA_FILE $ORIGINAL_BG_MARKOV_MODEL_FILE
#    echo -e "\nMarkov model generation for sorted original background file completed...\n"   
#
#    ### create zero order fasta markov model as a background file (useful for spamo analysis, if motif spacing analysis performed on original file rather than centered ones):
#    $MEME_SUITE_PATH/fasta-get-markov -m 0 $ORIGINAL_BG_FASTA_FILE $ORIGINAL_BG_ZERO_MARKOV_MODEL_FILE
#    echo -e "\nZero order Markov model generation for sorted original background file completed...\n"   
#
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
#   ### create zero order fasta markov model as a background file, useful esp. for spamo:
#   $MEME_SUITE_PATH/fasta-get-markov -m 0 $BG_FASTA_FILE $BG_ZERO_MARKOV_MODEL_FILE
#   echo -e "\nZero order Markov model generation for background file completed...\n"   

    ### Tomtom similarity analysis with 2 databases: 
#    $MEME_SUITE_PATH/tomtom -verbosity 1 -oc $TOMTOM_DIR -png -eps -min-overlap 4 -dist pearson -evalue -thresh 1 -no-ssc -bfile $BG_MARKOV_MODEL_FILE $MEME_DBFILE_PATH/meme.xml $MOTIF_DB_PATH/CIS-BP/Homo_sapiens.meme $MOTIF_DB_PATH/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme

#    $MEME_SUITE_PATH/tomtom -verbosity 1 -oc $TOMTOM_CISBP_DIR -png -eps -min-overlap 4 -dist pearson -evalue -thresh 1 -no-ssc -bfile $BG_MARKOV_MODEL_FILE $MEME_DBFILE_PATH/meme.xml $MOTIF_DB_PATH/CIS-BP/Homo_sapiens.meme 

    $MEME_SUITE_PATH/tomtom -verbosity 1 -oc $TOMTOM_COMP_CISBP_DIR -png -eps -min-overlap 4 -dist pearson -evalue -thresh 100 -no-ssc -bfile $BG_MARKOV_MODEL_FILE $MEME_DBFILE_PATH/meme.xml $MOTIF_DB_PATH/Homosapiens_custom_cisbp_allmotifs.meme

    $MEME_SUITE_PATH/tomtom -verbosity 1 -oc $TOMTOM_JASPAR2018_DIR -png -eps -min-overlap 4 -dist pearson -evalue -thresh 100 -no-ssc -bfile $BG_MARKOV_MODEL_FILE $MEME_DBFILE_PATH/meme.xml /gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme

    $MEME_SUITE_PATH/tomtom -verbosity 1 -oc $TOMTOM_JASPAR2016_DIR -png -eps -min-overlap 4 -dist pearson -evalue -thresh 100 -no-ssc -bfile $BG_MARKOV_MODEL_FILE $MEME_DBFILE_PATH/meme.xml /gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/JASPAR/JASPAR_CORE_2016_vertebrates.meme

    echo -e "\nTOMTOM analysis completed...\n"
    
#    $MEME_SUITE_PATH/fimo --max-stored-scores 150000 --parse-genomic-coord --verbosity 1 --oc $FIMO_DIR --bgfile $ORIGINAL_BG_MARKOV_MODEL_FILE $MEME_DBFILE_PATH/meme.xml  ${ORIGINAL_FASTA_FILE} 
#    echo -e "\nMEME-based FIMO analysis completed...\n"
#
#    $MEME_SUITE_PATH/centrimo --seqlen 500 --verbosity 1 --oc $CENTRIMO_DIR --neg $BG_FASTA_FILE --disc --bfile $BG_MARKOV_MODEL_FILE --local --score 5 --ethresh 10 $CENTERED_FASTA_FILE $MEME_DBFILE_PATH/meme.xml 
#    
#    echo -e "\nMEME-based CENTRIMO analysis completed...\n"
#
#    ### Make sure that you scan the motif within the original sequence and not on the 500bp centered fasta file, which is useful for centrimo and spamo purpose only:
#
#    # Both cisbp-jaspar together with meme motifs:
#    if (egrep -i -q "\b${TF_NAME_IDENTIFIER}\b" $MOTIF_DB_PATH/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme) && (egrep -i -q "\b${TF_NAME_IDENTIFIER}\b" $MOTIF_DB_PATH/CIS-BP/Homo_sapiens.meme); then 
#        echo -e "\n${TF_NAME_IDENTIFIER} : Present in both CISBP-JASPAR database\n"
#        echo -e "\nAnalyzing : $MEME_DBFILE_PATH meme file of ${TF_NAME_IDENTIFIER}\n"
#        # extract meme format for both cisbp and jaspar:
#        $MEME_SUITE_PATH/meme-get-motif -id $(egrep -i "\b${TF_NAME_IDENTIFIER}\b"  ${MOTIF_DB_PATH}/CIS-BP/Homo_sapiens.meme | head -1 | cut -f2 -d " ") ${MOTIF_DB_PATH}/CIS-BP/Homo_sapiens.meme > ${CISBP_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme
#
#        $MEME_SUITE_PATH/meme-get-motif -id $(egrep -i "\b${TF_NAME_IDENTIFIER}\b"  ${MOTIF_DB_PATH}/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme | head -1| cut -f2 -d " ") ${MOTIF_DB_PATH}/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme > ${JASPAR_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme
#
#        $MEME_SUITE_PATH/centrimo --seqlen 500 --verbosity 1 --oc $CENTRIMO_CISBP_JASPAR_DIR --neg $BG_FASTA_FILE --disc --bfile $BG_MARKOV_MODEL_FILE --local --score 5 --ethresh 10 $CENTERED_FASTA_FILE $MEME_DBFILE_PATH/meme.xml ${CISBP_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme ${JASPAR_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme 
#        echo -e "\nMEME-CISBP-JASPAR based Centrimo analysis completed...\n"
#    fi   
#
#    ### Centrimo analysis with equal-size ChIP-ed regions with 250bp up-downstream:
#    if egrep -i -q "\b${TF_NAME_IDENTIFIER}\b" $MOTIF_DB_PATH/CIS-BP/Homo_sapiens.meme; then 
#        $MEME_SUITE_PATH/meme-get-motif -id $(egrep -i "\b${TF_NAME_IDENTIFIER}\b"  ${MOTIF_DB_PATH}/CIS-BP/Homo_sapiens.meme | head -1 | cut -f2 -d " ") ${MOTIF_DB_PATH}/CIS-BP/Homo_sapiens.meme > ${CISBP_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme
#
#        $MEME_SUITE_PATH/fimo --max-stored-scores 150000 --parse-genomic-coord --verbosity 1 --oc $FIMO_CISBP_DIR --bgfile $ORIGINAL_BG_MARKOV_MODEL_FILE ${CISBP_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme ${ORIGINAL_FASTA_FILE} 
#
#        echo -e "\nCISBP based FIMO analysis completed...\n"
#        echo -e "\nAnalyzing : $MEME_DBFILE_PATH meme file of ${TF_NAME_IDENTIFIER}\n"
#        $MEME_SUITE_PATH/centrimo --seqlen 500 --verbosity 1 --oc $CENTRIMO_CISBP_DIR --neg $BG_FASTA_FILE --disc --bfile $BG_MARKOV_MODEL_FILE --local --score 5 --ethresh 10 $CENTERED_FASTA_FILE $MEME_DBFILE_PATH/meme.xml ${CISBP_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme 
#        echo -e "\nCISBP based Centrimo analysis completed...\n"
#    fi   
#
#    # For JASPAR:
#    if egrep -i -q "\b${TF_NAME_IDENTIFIER}\b" $MOTIF_DB_PATH/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme; then 
#        $MEME_SUITE_PATH/meme-get-motif -id $(egrep -i "\b${TF_NAME_IDENTIFIER}\b"  ${MOTIF_DB_PATH}/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme | head -1 | cut -f2 -d " ") ${MOTIF_DB_PATH}/JASPAR/JASPAR2018_CORE_vertebrates_non-redundant.meme > ${JASPAR_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme
#
#        $MEME_SUITE_PATH/fimo --max-stored-scores 150000 --parse-genomic-coord --verbosity 1 --oc $FIMO_JASPAR_DIR --bgfile $ORIGINAL_BG_MARKOV_MODEL_FILE ${JASPAR_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme ${ORIGINAL_FASTA_FILE} 
#        echo -e "\nJASPAR based FIMO analysis completed...\n"
#        echo -e "\nAnalyzing : $MEME_DBFILE_PATH meme file of ${TF_NAME_IDENTIFIER}\n"
#
#        $MEME_SUITE_PATH/centrimo --seqlen 500 --verbosity 1 --oc $CENTRIMO_JASPAR_DIR --neg $BG_FASTA_FILE --disc --bfile $BG_MARKOV_MODEL_FILE --local --score 5 --ethresh 10 $CENTERED_FASTA_FILE $MEME_DBFILE_PATH/meme.xml ${JASPAR_MEME_DIR}/${TF_NAME_IDENTIFIER}.meme 
#        echo -e "\nJASPAR based Centrimo analysis completed...\n"
#    fi   
#
fi

echo -e "\nMotif analysis completed!!!\n"




