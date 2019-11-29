#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing
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

if [ -z "$OUTPUT_DIR" ]; then OUTPUT_DIR=$(get_output_dir $LIBRARY);  fi
if [ ! -d "$OUTPUT_DIR" ]; then 
    mkdir -p $OUTPUT_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $OUTPUT_DIR. Exiting"; exit 1; fi
fi

if [ -z "$TEMP_DIR" ]; then TEMP_DIR=$(get_temp_dir); fi
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p $TEMP_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $TEMP_DIR. Exiting"; exit 1; fi
fi

### This script will created pooled replicates, pooled pseudoreplicates, and pooled controls ###
### It will need Rep1 SL#, Rep2 SL#, Control1 SL#, and Control2 SL# passed in on the command ###
### line.                                                                                    ###
### This script assumes FastqToTagAlign has been run on the above 4 samples.                 ###

if [ $# -ne 5 ]; then echo "Usage $0 <Rep1 SL#> <Rep2 SL#> <Control1 SL#> <Control2 SL#> <Output Dir>"; exit 1; fi

REP1=$1
REP2=$2
CONTROL1=$3
CONTROL2=$4
MY_OUTPUT_DIR=$5

### TEMP FOR TESTING ###
REP1_DIR=$BASE_DIR/Analysis/$REP1
REP2_DIR=$BASE_DIR/Analysis/$REP2
CONTROL1_DIR=$BASE_DIR/Analysis/$CONTROL1
CONTROL2_DIR=$BASE_DIR/Analysis/$CONTROL2
#OUTPUT_DIR=$BASE_DIR/IDR_${REP1}_${REP2}
OUTPUT_DIR=$MY_OUTPUT_DIR
mkdir -p $OUTPUT_DIR
### END TEMP FOR TESTING ###

# Get inputs together #
REP1_TA_FILE=$REP1_DIR/$REP1.filt.nodup.srt.SE.tagAlign.gz
REP2_TA_FILE=$REP2_DIR/$REP2.filt.nodup.srt.SE.tagAlign.gz

CONTROL1_TA_FILE=$CONTROL1_DIR/$CONTROL1.filt.nodup.srt.SE.tagAlign.gz
CONTROL2_TA_FILE=$CONTROL2_DIR/$CONTROL2.filt.nodup.srt.SE.tagAlign.gz

REP1_PR1_TA_FILE=$REP1_DIR/$REP1.filt.nodup.SE.pr1.tagAlign.gz
REP1_PR2_TA_FILE=$REP1_DIR/$REP1.filt.nodup.SE.pr2.tagAlign.gz

REP2_PR1_TA_FILE=$REP2_DIR/$REP2.filt.nodup.SE.pr1.tagAlign.gz
REP2_PR2_TA_FILE=$REP2_DIR/$REP2.filt.nodup.SE.pr2.tagAlign.gz

# Name output files #
DATASET_PREFIX=${REP1}_${REP2}
POOLED_TA_FILE=$OUTPUT_DIR/$DATASET_PREFIX.Rep0.tagAlign.gz
POOLED_CONTROL_FILE=$OUTPUT_DIR/$DATASET_PREFIX.Rep0.control.tagAlign.gz
PPR1_TA_FILE=$OUTPUT_DIR/$DATASET_PREFIX.Rep0.pr1.tagAlign.gz
PPR2_TA_FILE=$OUTPUT_DIR/$DATASET_PREFIX.Rep0.pr2.tagAlign.gz

CMD="zcat $REP1_TA_FILE $REP2_TA_FILE | gzip -c > $POOLED_TA_FILE"
run_cmd "$CMD" "$POOLED_TA_FILE"

if [ "$CONTROL1_TA_FILE" = "$CONTROL2_TA_FILE" ]; then
    CMD="cp $CONTROL1_TA_FILE $POOLED_CONTROL_FILE"
else
    CMD="zcat $CONTROL1_TA_FILE $CONTROL2_TA_FILE | gzip -c > $POOLED_CONTROL_FILE"
fi
run_cmd "$CMD" "$POOLED_CONTROL_FILE"

CMD="zcat $REP1_PR1_TA_FILE $REP2_PR1_TA_FILE | gzip -c > $PPR1_TA_FILE"
run_cmd "$CMD" "$PPR1_TA_FILE"

CMD="zcat $REP1_PR2_TA_FILE $REP2_PR2_TA_FILE | gzip -c > $PPR2_TA_FILE"
run_cmd "$CMD" "$PPR2_TA_FILE"

log_msg info "Pooling finished..."
