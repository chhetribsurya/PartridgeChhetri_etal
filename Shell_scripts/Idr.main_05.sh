#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing
#OUTPUT_DIR=$BASE_DIR/Analysis/$LIBRARY
#TEMP_DIR=$BASE_DIR/TEMP
#REFERENCE_DIR=$BASE_DIR/../reference

#source /gpfs/gpfs1/home/schhetri/python/anaconda_python_version.sh
#source activate python3

#export PATH="/gpfs/gpfs1/software/python3.3/bin/:$PATH"
#SOURCE_DIR="/opt/HAIB/myerslab/etc"
SOURCE_DIR="/gpfs/gpfs2/software/HAIB/myerslab/etc"
export PATH=/gpfs/gpfs2/software/python-3.6.0/bin:${PATH}
PYTHONPATH="" #unset PYTHONPATH
#source /gpfs/gpfs1/home/schhetri/python/anaconda_python_version.sh
#source /gpfs/gpfs1/home/schhetri/Tools/python3-venv/bin/activate

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

### This script will be passed a pair of narrow peak files to run IDR against. It will also need the ###
### pooled peak file and a blacklist. Finally, it will need a value to pass to --rank                ###

if [ $# -ne 5 ]; then echo "Usage: $0 <REP1 peaks file> <REP2 peaks file> <Pooled Peaks file> <Rank param> <Output dir>"; exit 1; fi

REP1_PEAK_FILE=$1
REP2_PEAK_FILE=$2
REP1=$(basename $REP1_PEAK_FILE .narrowPeak.gz)
REP2=$(basename $REP2_PEAK_FILE .narrowPeak.gz)
REP1_VS_REP2="${REP1}_VS_${REP2}"
POOLED_PEAK_FILE=$3
RANK=$4
OUTPUT_DIR=$5
IDR_OUTPUT=$OUTPUT_DIR/${REP1_VS_REP2}.IDR0.narrowPeak.gz
IDR_THRESHOLD_FILE=$OUTPUT_DIR/${REP1_VS_REP2}.IDR0.05.narrowPeak.gz
IDR_FILT_THRESHOLD_FILE=$OUTPUT_DIR/${REP1_VS_REP2}.IDR0.05.filt.narrowPeak.gz

BLACKLIST=$REFERENCE_DIR/wgEncodeDacMapabilityConsensusExcludable.bed.gz
IDR_THRESH=0.05

#IDR_CMD=/gpfs/gpfs2/sdi/software/idr-2.0.2/bin/idr
IDR_CMD=/gpfs/gpfs1/home/schhetri/Tools/idr-2.0.2/bin/idr
IDR_PEAK_DIR=$OUTPUT_DIR/idr_peaks

CMD="$IDR_CMD --samples ${REP1_PEAK_FILE} ${REP2_PEAK_FILE} --peak-list ${POOLED_PEAK_FILE} --input-file-type narrowPeak --output-file ${IDR_OUTPUT} --rank $RANK --soft-idr-threshold ${IDR_THRESH} --plot"
run_cmd "$CMD" "$IDR_OUTPUT"

IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')

CMD="awk 'BEGIN{OFS=\"\t\"} \$12>='\"${IDR_THRESH_TRANSFORMED}\"' {print \$1,\$2,\$3,\$4,\$5,\$6,\$7,\$8,\$9,\$10}' ${IDR_OUTPUT} | sort | uniq | egrep -v '^chrM\b' | sort -k7n,7n | gzip -c > ${IDR_THRESHOLD_FILE}"
run_cmd "$CMD" "$IDR_THRESHOLD_FILE"

NPEAKS_IDR=$(zcat $IDR_THRESHOLD_FILE | wc -l)
log_msg info "Number of peaks passing IDR: $NPEAKS_IDR"

BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin/bedtools"
zcat $IDR_THRESHOLD_FILE | ${BEDTOOLS_PATH} intersect -v -a - -b ${BLACKLIST} | gzip -c > $IDR_FILT_THRESHOLD_FILE

#BEDTOOLS=$(get_software_dir bedtools2-2.20.0/bin/bedtools)
#CMD="zcat $IDR_THRESHOLD_FILE | $BEDTOOLS intersect -v -a - -b $BLACKLIST | gzip -c > $IDR_FILT_THRESHOLD_FILE"
#run_cmd "$CMD" "$IDR_FILT_THRESHOLD_FILE"
#zcat $IDR_THRESHOLD_FILE | ${BEDTOOLS_PATH} intersect -v -a - -b <(zcat -f ${BLACKLIST} | gzip -c > $IDR_FILT_THRESHOLD_FILE

if [[ ! -d $IDR_PEAK_DIR ]]; then
    mkdir -p $IDR_PEAK_DIR
fi

if [[ -f $IDR_FILT_THRESHOLD_FILE ]]; then 
    cp $IDR_FILT_THRESHOLD_FILE $IDR_PEAK_DIR
fi

if [[ -f $IDR_THRESHOLD_FILE ]]; then
    cp $IDR_THRESHOLD_FILE $IDR_PEAK_DIR
fi

if [[ -f $IDR_OUTPUT ]]; then
    cp $IDR_OUTPUT $IDR_PEAK_DIR
fi

if [[ -f $OUTPUT_DIR/idr*truerep.out ]]; then
    cp $OUTPUT_DIR/idr*truerep.out $IDR_PEAK_DIR 
fi

log_msg info "IDR Complete."
