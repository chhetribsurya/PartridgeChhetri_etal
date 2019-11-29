#!/bin/bash

export RUN_PATH=`pwd`
#To overcome the issue of exporting the bash array variable, so directly: (Read lines in file into an bash array)
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


for i in "${!LIB_REP1[@]}"; do

        R1=${LIB_REP1[$i]}
        R2=${LIB_REP2[$i]}
	IDR_DIR=$BASE_DIR/IDR_${R1}_${R2}
	echo -e "\ndir for logfile : $IDR_DIR/motif_idr.${R1}_${R2}.truerep.out" 
	#echo -e "\ncheck all idr peaks agnostic of cutoff: $IDR_DIR/${R1}.filt.nodup.srt.SE_VS_${R2}.filt.nodup.srt.SE.IDR0.narrowPeak.gz\n"
	#echo -e "\ncheck idr passed peaks with cutoff 0.02 (Real peaks for the analysis) : $IDR_DIR/${R1}.filt.nodup.srt.SE_VS_${R2}.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz\n"
	
	#rm $IDR_DIR/motif_corrected_idr.${R1}_${R2}.truerep.out
	# Motif analysis using MeMe
	#bsub $BSUB_OPTIONS -o "$IDR_DIR/motif_corrected_idr.${R1}_${R2}.truerep.out" -J "Motif analysis from IDR of $R1 $R2 " $RUN_PATH/meme_peaks_motif_analysis.sh $IDR_DIR/${R1}.filt.nodup.srt.SE_VS_${R2}.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz $R1 $R2
	bsub $BSUB_OPTIONS -o "./motif_corrected_idr.${R1}_${R2}.truerep.out" -J "Motif analysis from IDR of $R1 $R2 " $RUN_PATH/meme_peaks_motif_analysis.sh $IDR_DIR/${R1}.filt.nodup.srt.SE_VS_${R2}.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz $R1 $R2

done

