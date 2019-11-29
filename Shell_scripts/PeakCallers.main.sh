#!/bin/bash

#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing
#OUTPUT_DIR=$BASE_DIR/Analysis/$LIBRARY
#TEMP_DIR=$BASE_DIR/TEMP
#REFERENCE_DIR=$BASE_DIR/../reference

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
if [ -z "$LSB_JOBID" ]; then echo "Please run on a compute node. Exiting"; exit 1; fi

if [ -z "$TEMP_DIR" ]; then TEMP_DIR=$(get_temp_dir); fi
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p $TEMP_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $TEMP_DIR. Exiting"; exit 1; fi
fi

### This script will be passed a ChIP tagalign file and a Control tagAlign file on the command line ###
### The reason for this is we have to run Peak Calling against a ton of different replicates and    ###
###   pseudoreplicates.                                                                             ###
### For this reason, the output directory must also be passed in. We'll basename that param for the ***
###   file prefix. Also we'll need a fragment length                                                      ###

if [ $# -ne 4 ]; then echo "Usage: $0 <CHIP TA FILE> <CONTROL TA FILE> <OUTPUT DIR> <FRAG LEN>"; exit 1; fi

CHIP_TA_FILE=$1
CONTROL_TA_FILE=$2
OUTPUT_DIRECTORY=$3
FRAGLEN=$4

CHIP_TA_PREFIX=$(basename $CHIP_TA_FILE .tagAlign.gz)

### For now just do SPP; will change this around to get other peak callers up and running ###
SPP="Y"
GEM="N"
PEAKSEQ="N"
### END For now ...                                                                       ###

#export R_LIBS_SITE=/gpfs/gpfs1/software/R-site-packages-3.0
export R_LIBS_SITE=/gpfs/gpfs2/software/c7R-libs

# Give the option to turn off peak callers if desired
if [ -z "$SPP" ]; then SPP="Y"; fi
if [ -z "$GEM" ]; then GEM="Y"; fi
if [ -z "$PEAKSEQ" ]; then PEAKSEQ="Y"; fi

if [ "$SPP" != "N" ]; then
    SPP_NODUPS=$(get_software_dir spp-1.1/run_spp_nodups.R)
    CMD="/gpfs/gpfs2/software/R-3.3.3/bin/Rscript $SPP_NODUPS -c=$CHIP_TA_FILE -i=$CONTROL_TA_FILE -npeak=300000 -odir=$OUTPUT_DIRECTORY -speak=$FRAGLEN -savr -savn=$OUTPUT_DIRECTORY/$CHIP_TA_PREFIX.narrowPeak -savp -rf -out=$OUTPUT_DIRECTORY/$CHIP_TA_PREFIX.ccscores"
    run_cmd "$CMD" "$OUTPUT_DIRECTORY/$CHIP_TA_PREFIX.narrowPeak.gz $OUTPUT_DIRECTORY/$CHIP_TA_PREFIX.tagAlign.pdf $OUTPUT_DIRECTORY/$CHIP_TA_PREFIX.ccscores"
fi

if [ "$GEM" != "N" ]; then
    GEM_JAR=$(get_software_dir gem/gem.jar)
    CMD="zcat $CHIP_TA_FILE > $TEMP_DIR/$CHIP_TA_PREFIX.tagAlign"
    run_cmd "$CMD" "$TEMP_DIR/$CHIP_TA_PREFIX.tagAlign"

    #CMD="java -Xmx15G -jar $GEM_JAR --g $REFERENCE_DIR/hg19-male.gem.info --d $(dirname $GEM_JAR)/Read_Distribution_default.txt --s 2400000000 --expt $TEMP_DIR/$CHIP_TA_PREFIX.tagAlign --ctrl $CONTROL_TA_FILE --f BED --out $OUTPUT_DIRECTORY/$CHIP_TA_PREFIX --genome 
fi

if [ "$PEAKSEQ" != "N" ]; then
	echo "Note: only SPP called chip peaks processsed but Peakseq not called"

fi
