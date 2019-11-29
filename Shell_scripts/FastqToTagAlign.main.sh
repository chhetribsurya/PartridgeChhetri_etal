#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing
#OUTPUT_DIR=$BASE_DIR/Analysis/$LIBRARY
#TEMP_DIR=$BASE_DIR/TEMP
#REFERENCE_DIR=$BASE_DIR/../reference

SOURCE_DIR="/gpfs/gpfs2/software/HAIB/myerslab/etc"
#SOURCE_DIR="/opt/HAIB/myerslab/etc"

### Mandatory sourcing of bashrc for necessary environment variables. ###
if [ -e $SOURCE_DIR/bashrc ]; then 
    . $SOURCE_DIR/bashrc
else echo "[fatal] - Could not find myerslab bashrc file. Exiting"; exit 1; fi

### Mandatory sourcing of functions to get helper functions (like call_cmd). ###
if [ -e $SOURCE_DIR/functions ]; then
    . $SOURCE_DIR/functions
else echo "[fatal] - Could not find functions file. Exiting"; exit 1; fi

#export R_LIBS_SITE=/gpfs/gpfs1/software/R-site-packages
export R_LIBS_SITE=/gpfs/gpfs2/software/c7R-libs

### Verify we are not running on the head node. ###
if [ -z "$LSB_JOBID" ]; then echo "Please run on a compute node. Exiting"; exit 1; fi

if [ -z "$LIBRARY" ];            then empty_param_quit "LIBRARY"; fi
if [ -z "$GENOME" ];             then empty_param_quit "GENOME"; fi

### Verify library directory exists. ###
if [ -z "$LIBRARY_DIR" ]; then LIBRARY_DIR=$(get_library_dir $LIBRARY); fi
if [ ! -d "$LIBRARY_DIR" ]; then echo "Library directory does not exist: $LIBRARY_DIR. Exiting"; exit 1; fi

if [ -z "$OUTPUT_DIR" ]; then OUTPUT_DIR=$(get_output_dir $LIBRARY);  fi
if [ ! -d "$OUTPUT_DIR" ]; then 
    mkdir -p $OUTPUT_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $OUTPUT_DIR. Exiting"; exit 1; fi
fi

### Set up temp directory variable ###
if [ -z "$TEMP_DIR" ]; then TEMP_DIR=$(get_temp_dir); fi
if [ ! -d "$TEMP_DIR" ]; then
    mkdir -p $TEMP_DIR
    if [ $? -ne "0" ]; then echo "Could not create output dir: $TEMP_DIR. Exiting"; exit 1; fi
fi


export LOGFILE_NAME=$OUTPUT_DIR/$LIBRARY.log

#NTHREADS=8
NTHREADS=6

log_msg info "Beginning Fastq to TagAlign"
log_msg info "    LIBRARY:              $LIBRARY"
log_msg info "    GENOME:               $GENOME"
log_msg info "    LIBRARY_DIR:          $LIBRARY_DIR"
log_msg info "    OUTPUT_DIR:           $OUTPUT_DIR"

log_msg info "Consolidating fastq files..."
# Set up output files
#FASTQ_FILE_1=$TEMP_DIR/$LIBRARY.fastq.gz
#SAI_FILE_1=$TEMP_DIR/$LIBRARY.sai
#RAW_BAM_PREFIX=$LIBRARY.raw.srt
#RAW_BAM_FILE=$OUTPUT_DIR/$RAW_BAM_PREFIX.bam
#RAW_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$RAW_BAM_PREFIX.flagstat.qc

# Set up output files
FASTQ_FILE_1=$TEMP_DIR/$LIBRARY.fastq.gz
SAI_FILE_1=$TEMP_DIR/$LIBRARY.sai
SAM_FILE=$TEMP_DIR/$LIBRARY.raw.sam
BAM_FILE=$TEMP_DIR/$LIBRARY.raw.bam
RAW_BAM_PREFIX=$LIBRARY.raw.srt
RAW_BAM_FILE=$OUTPUT_DIR/$RAW_BAM_PREFIX.bam
RAW_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$RAW_BAM_PREFIX.flagstat.qc

CMD="cat $LIBRARY_DIR/*fastq.gz > $FASTQ_FILE_1"
run_cmd "$CMD" "$FASTQ_FILE_1"

BWA_SOFTWARE=$(get_software_dir bwa-0.7.12/bwa)

log_msg info "Aligning with bwa aln..."
CMD="$BWA_SOFTWARE aln -q 5 -l 32 -k 2 -t $NTHREADS $REFERENCE_DIR/$GENOME.fa $FASTQ_FILE_1 > $SAI_FILE_1"
run_cmd "$CMD" "$SAI_FILE_1"

#CMD="$BWA_SOFTWARE samse $REFERENCE_DIR/$GENOME.fa $SAI_FILE_1 $FASTQ_FILE_1 | samtools view -Su - | samtools sort - $OUTPUT_DIR/$RAW_BAM_PREFIX"
$BWA_SOFTWARE samse -f ${SAM_FILE} $REFERENCE_DIR/$GENOME.fa $SAI_FILE_1 $FASTQ_FILE_1 
samtools view -bS ${SAM_FILE} > ${BAM_FILE}
samtools sort ${BAM_FILE} -o $RAW_BAM_FILE
##rm $BAM_FILE
##rm $SAI_FILE_1

CMD="samtools flagstat $RAW_BAM_FILE > $RAW_BAM_FILE_MAPSTATS"
run_cmd "$CMD" "$RAW_BAM_FILE_MAPSTATS"

log_msg info "Post-alignment filtering step..."
# Set up filter bam output files
FILT_BAM_PREFIX=$LIBRARY.filt.srt
FILT_BAM_FILE=$OUTPUT_DIR/$FILT_BAM_PREFIX.bam
DUP_FILE_QC=$OUTPUT_DIR/$FILT_BAM_PREFIX.dup.qc

MAPQ_THRESH=30
CMD="samtools view -F 1804 -q $MAPQ_THRESH -b $RAW_BAM_FILE > $FILT_BAM_FILE"
run_cmd "$CMD" "$FILE_BAM_FILE"

TMP_FILT_BAM_FILE=$TEMP_DIR/$FILT_BAM_PREFIX.dupmark.bam
PICARD_JAR=$(get_software_dir picard-tools-1.88/MarkDuplicates.jar)
CMD="java -Xmx4G -jar $PICARD_JAR INPUT=$FILT_BAM_FILE OUTPUT=$TMP_FILT_BAM_FILE METRICS_FILE=$DUP_FILE_QC VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false"
run_cmd "$CMD" "$TMP_FILT_BAM_FILE $DUP_FILE_QC"
mv $TMP_FILT_BAM_FILE $FILT_BAM_FILE

# Set up final bam output files
FINAL_BAM_PREFIX=$LIBRARY.filt.nodup.srt
FINAL_BAM_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.bam
FINAL_BAM_INDEX_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.bai
FINAL_BAM_FILE_MAPSTATS=$OUTPUT_DIR/$FINAL_BAM_PREFIX.flagstat.qc

CMD="samtools view -F 1804 -b $FILT_BAM_FILE > $FINAL_BAM_FILE"
run_cmd "$CMD" "$FINAL_BAM_FILE"
CMD="samtools index $FINAL_BAM_FILE $FINAL_BAM_INDEX"
run_cmd "$CMD" "$FINAL_BAM_INDEX"
CMD="samtools flagstat $FINAL_BAM_FILE > $FINAL_BAM_FILE_MAPSTATS"
run_cmd "$CMD" "$FINAL_BAM_FILE_MAPSTATS"

# Set up library complexity output file
PBC_FILE_QC=$OUTPUT_DIR/$FINAL_BAM_PREFIX.pbc.qc

BEDTOOLS=$(get_software_dir bedtools2-2.20.0/bin/bedtools)
CMD="$BEDTOOLS bamtobed -i $FILT_BAM_FILE | awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$6}' | grep -v 'chrM' | sort | uniq -c | awk 'BEGIN{mt=0;m0=0;m1=0;m2=0} (\$1==1){m1=m1+1} (\$1==2){m2=m2+1} {m0=m0+1} {mt=mt+1} END{printf \"%d\t%d\t%d\t%d\t%f\t%f\t%f\n\",mt,m0,m1,m2,m0/mt,m1/m0,m1/m2}' > $PBC_FILE_QC"
run_cmd "$CMD" "$PBC_FILE_QC"

rm $FILT_BAM_FILE

log_msg info "Create tag align files..."
# Set up tag align output files
NREADS=15000000
FINAL_TA_FILE=$OUTPUT_DIR/$FINAL_BAM_PREFIX.SE.tagAlign.gz
SUBSAMPLED_TA_FILE=$OUTPUT_DIR/$LIBRARY.filt.nodup.sample.$(( NREADS / 1000000 )).SE.tagAlign.gz

CMD="$BEDTOOLS bamtobed -i $FINAL_BAM_FILE | awk 'BEGIN{OFS=\"\t\"}{\$4=\"N\";\$4=\"1000\";print \$0}' | gzip -c > $FINAL_TA_FILE"
run_cmd "$CMD" "$FINAL_TA_FILE"

CMD="zcat $FINAL_TA_FILE | grep -v 'chrM' | shuf -n $NREADS | gzip -c > $SUBSAMPLED_TA_FILE"
run_cmd "$CMD" "$SUBSAMPLED_TA_FILE"

log_msg info "Calcuate cross-correlation QC scores..."
# Set up CC QC output files
CC_SCORES_FILE=$SUBSAMPLED_TA_FILE.cc.qc
CC_PLOT_FILE=$SUBSAMPLED_TA_FILE.cc.plot.pdf


### There seems to be bug on this line of code, since it outputs the same $CC_SCORES_FILE if re-run, and that would lead to repetition of the
### same line, and thus when you cut that line with $(cut -f 3 $CC_SCORES_FILE) then you are passing the 2 fragment length. This would kill
### the program which can only acccomodate only single fragment length argument. So, maybe you could just integrate $(cut -f 3 $CC_SCORES_FILE | tr "\n" "\t" | cut -f 1)
### in submit_peak_calls.sh script ; this would avoid the duplication of the same fragment. 
#SPP_NODUPS=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/spp_install/phantompeakqualtools/run_spp_nodups.R
SPP_NODUPS=/gpfs/gpfs1/home/schhetri/spp_install/phantompeakqualtools/run_spp_nodups.R
CMD="Rscript $SPP_NODUPS -c=$SUBSAMPLED_TA_FILE -p=$NTHREADS -rf -filtchr=chrM -savp=$CC_PLOT_FILE -out=$CC_SCORES_FILE"
run_cmd "$CMD" "$CC_SCORES_FILE $CC_PLOT_FILE"

sed -ri 's/,[^\t]+//g' $CC_SCORES_FILE

log_msg info "Creating self-psuedoreplicates..."
# Set up PR output files
PR_PREFIX=$LIBRARY.filt.nodup
PR1_TA_FILE=$OUTPUT_DIR/$PR_PREFIX.SE.pr1.tagAlign.gz
PR2_TA_FILE=$OUTPUT_DIR/$PR_PREFIX.SE.pr2.tagAlign.gz

nlines=$( zcat $FINAL_TA_FILE | wc -l )
nlines=$(( (nlines + 1) / 2 ))

CMD="zcat $FINAL_TA_FILE | shuf | split -d -l $nlines - $TEMP_DIR/$PR_PREFIX"
run_cmd "$CMD" "$TEMP_DIR/${PR_PREFIX}00 $TEMP_DIR/${PR_PREFIX}01"

CMD="gzip -c $TEMP_DIR/${PR_PREFIX}00 > $PR1_TA_FILE"
run_cmd "$CMD" "$PR1_TA_FILE"
CMD="gzip -c $TEMP_DIR/${PR_PREFIX}01 > $PR2_TA_FILE"
run_cmd "$CMD" "$PR2_TA_FILE"

log_msg info "Tag align script complete. Output is in $OUTPUT_DIR"

echo ""
echo "starting the script for the pooling of tag aligns:"
### Calling the script for pooling of the true tag align and control tag align:
#bsub -We 24:00 -n 1 -R span[hosts=1] -J "Call Pooling data" -o $LOG_DIR/pool_data.out $RUN_PATH/submit_PoolDataSets.main.sh
