#!/bin/bash

export RUN_PATH=`pwd`

#if [ $# -ne 2 ]; then echo "Usage: $0 <LIBRARY> <GENOME>"; exit 1; fi

#LIBRARY=$1
#GENOME=$2

echo "starting the script"
echo "$GENOME"
echo "$BASE_DIR"

#Since bash array variable cannot be exported so the following array variable is not imported:
echo "library ${LIB_REP1[@]}" #Thus this line won't be executed 

#readarray ALL_LIB < ./all_lib.txt
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
#sort -u batch_VI_all_tf.txt | grep -v "SL154817"

#Reading or geting unique values from an bash array, and then note it's reconverted into bash array by placing array bracket ( ):
#UNIQ_ALL_LIB=( $(echo ${ALL_LIB[@]} | tr " " "\n" | sort -u | tr "\n" " ") )


#for i in {1..40}; do
#	echo $i
#	sleep 1
#done

for i in "${!UNIQ_ALL_LIB[@]}";do 
	export LIBRARY=${UNIQ_ALL_LIB[$i]} 
	export LIBRARY_DIR=$LIB_PATH/$LIBRARY
	export OUTPUT_DIR=$BASE_DIR/Analysis/$LIBRARY
	export GENOME="hg19-male"
	echo "processing $LIBRARY_DIR"
	
	if [[ ! -d $OUTPUT_DIR ]];then
		mkdir -p $OUTPUT_DIR
		if [[ $? -ne "0" ]]; then echo "Could not create the output dir: $OUTPUT_DIR. Exiting"; exit 1; fi
	fi

	bsub $BSUB_OPTIONS -J "step1a_${LIBRARY}_align" -o $LOG_DIR/${LIBRARY}_fastq_align.out $RUN_PATH/FastqToTagAlign.main.sh ${LIBRARY} ${GENOME}
done



