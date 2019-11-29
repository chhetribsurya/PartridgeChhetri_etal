#!/bin/bash

#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing
#BSUB_OPTIONS="-n 8 -R span[hosts=1]"

#if [ $# -ne 5 ]; then echo "Usage: $0 REP1 REP2 CONTROL1 CONTROL2 OUTPUT_DIR"; exit 1; fi

#REP1=$1
#REP2=$2
#CONTROL1=$3
#CONTROL2=$4
#MY_OUTPUT_DIR=$5
#MY_OUTPUT_DIR=$BASE_DIR/IDR_${REP1}_${REP2}/

#So, to overcome the issue of exporting the bash array variable, so directly: (Read lines in file into an bash array)

readarray LIB_REP1 < $REP1_FILE_NAME
readarray LIB_REP2 < $REP2_FILE_NAME
readarray LIB_CONTROL1 < $CONTROL1_FILE_NAME
readarray LIB_CONTROL2 < $CONTROL2_FILE_NAME
readarray ALL_LIB < $ALL_TF_FILE_NAME

#remove the space character from the array elements:
LIB_REP1=( $(echo ${LIB_REP1[@]} | tr " " "\n" | tr "\n" " "))
LIB_REP2=( $(echo ${LIB_REP2[@]} | tr " " "\n" | tr "\n" " "))
LIB_CONTROL1=( $(echo ${LIB_CONTROL1[@]} | tr " " "\n" | tr "\n" " "))
LIB_CONTROL2=( $(echo ${LIB_CONTROL2[@]} | tr " " "\n" | tr "\n" " "))
UNIQ_ALL_LIB=( $(echo ${ALL_LIB[@]} | tr " " "\n" | sort -u | tr "\n" " ") )

#echo ${LIB_REP1[@]}

for i in "${!LIB_REP1[@]}"; do

    REP1=${LIB_REP1[$i]}
    REP2=${LIB_REP2[$i]}
    CONTROL1=${LIB_CONTROL1[$i]}
    CONTROL2=${LIB_CONTROL2[$i]}
	MY_OUTPUT_DIR=$BASE_DIR/IDR_${REP1}_${REP2}
	
	#CURRENT_REP1=${LIB_REP1[$i]}
	#CURRENT_REP2=${LIB_REP2[$i]}
	#CURRENT_CONTROL1=${LIB_CONTROL1[$i]}
	#CURRENT_CONTROL2=${LIB_CONTROL2[$i]}
	#LOG_NAME=$(echo ${CURRENT_REP1}_${CURRENT_REP2}_${CURRENT_CONTROL1}_${CURRENT_CONTROL2}_pool.out | tr " " "_")
	
	echo -e "\nStarting the script for pooling of chip TA align/true TA align and control TA align\n"
	$RUN_PATH/PoolDataSets.main.sh ${REP1} ${REP2} ${CONTROL1} ${CONTROL2} ${MY_OUTPUT_DIR}

	echo -e "\nNow starting the Peak caller script for the spp peak call\n"
	echo -e "Processing ${REP1}_${REP2}_${CONTROL1}_${CONTROL2} library combination\n"
	#echo $REP1.filt.nodup.srt.SE.tagAlign.gz

	FRAGLEN_REP1=$(cut -f 3 $BASE_DIR/Analysis/$REP1/$REP1.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc | tr "\n" "\t"| cut -f 1)
	#FRAGLEN_REP1=$(cut -f 3 $BASE_DIR/Analysis/$REP1/$REP1.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc)
	FRAGLEN_REP2=$(cut -f 3 $BASE_DIR/Analysis/$REP2/$REP2.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc | tr "\n" "\t"| cut -f 1)
	#FRAGLEN_REP2=$(cut -f 3 $BASE_DIR/Analysis/$REP2/$REP2.filt.nodup.sample.15.SE.tagAlign.gz.cc.qc)

	if [ $FRAGLEN_REP1 -gt $FRAGLEN_REP2 ]; then
	    FRAGLEN_POOL=$FRAGLEN_REP1
	else
	    FRAGLEN_POOL=$FRAGLEN_REP2
	fi

	# True Rep 1 vs control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP1}_true.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $BASE_DIR/Analysis/$REP1/$REP1.filt.nodup.srt.SE.tagAlign.gz $BASE_DIR/Analysis/$CONTROL1/$CONTROL1.filt.nodup.srt.SE.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_REP1

	# True Rep 2 vs control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP2}_true_test.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $BASE_DIR/Analysis/$REP2/$REP2.filt.nodup.srt.SE.tagAlign.gz $BASE_DIR/Analysis/$CONTROL2/$CONTROL2.filt.nodup.srt.SE.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_REP2

	# Rep 1 Self PR 1 vs control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP1}_pr1.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $BASE_DIR/Analysis/$REP1/$REP1.filt.nodup.SE.pr1.tagAlign.gz $BASE_DIR/Analysis/$CONTROL1/$CONTROL1.filt.nodup.srt.SE.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_REP1

	# Rep 1 Self PR 2 vs control
 	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP1}_pr2.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $BASE_DIR/Analysis/$REP1/$REP1.filt.nodup.SE.pr2.tagAlign.gz $BASE_DIR/Analysis/$CONTROL1/$CONTROL1.filt.nodup.srt.SE.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_REP1

	# Rep 2 Self PR 1 vs control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP2}_pr1_test.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $BASE_DIR/Analysis/$REP2/$REP2.filt.nodup.SE.pr1.tagAlign.gz $BASE_DIR/Analysis/$CONTROL2/$CONTROL2.filt.nodup.srt.SE.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_REP2

	# Rep 2 Self PR 2 vs control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP2}_pr2_test.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $BASE_DIR/Analysis/$REP2/$REP2.filt.nodup.SE.pr2.tagAlign.gz $BASE_DIR/Analysis/$CONTROL2/$CONTROL2.filt.nodup.srt.SE.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_REP2

	# Pooled True Rep vs Pooled control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP1}_${REP2}_pooled_true.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $MY_OUTPUT_DIR/${REP1}_${REP2}.Rep0.tagAlign.gz $MY_OUTPUT_DIR/${REP1}_${REP2}.Rep0.control.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_POOL

	# Pooled Psuedorep 1 vs pooled control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP1}_${REP2}_pooled_pr1.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $MY_OUTPUT_DIR/${REP1}_${REP2}.Rep0.pr1.tagAlign.gz $MY_OUTPUT_DIR/${REP1}_${REP2}.Rep0.control.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_POOL

	# Pooled Pseduorep 2 vs pooled control
	bsub $BSUB_OPTIONS -o "$MY_OUTPUT_DIR/${REP1}_${REP2}_pooled_pr2.out" -J "${REP1}_${REP2}_${CONTROL1}_${CONTROL2} Peak calling" ./PeakCallers.main.sh $MY_OUTPUT_DIR/${REP1}_${REP2}.Rep0.pr2.tagAlign.gz $MY_OUTPUT_DIR/${REP1}_${REP2}.Rep0.control.tagAlign.gz $MY_OUTPUT_DIR $FRAGLEN_POOL

done
