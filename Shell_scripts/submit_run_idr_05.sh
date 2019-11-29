#!/bin/bash

#source /gpfs/gpfs1/home/schhetri/python/anaconda_python_version3.sh
#source activate python3

#BASE_DIR=/gpfs/gpfs1/home/snewberry/interrupt/encode/spp_and_idr/idr_testing

#BSUB_OPTIONS="-n 8 -R span[hosts=1]"
#if [ $# -ne 3 ]; then echo "Usage $0: <Rep1 SL#> <Rep2 SL#> <IDR directory name>"; exit 1; fi

#R1=$1
#R2=$2
#IDR_DIR=$3

#So, to overcome the issue of exporting the bash array variable, so directly: (Read lines in file into an bash array)
#readarray LIB_REP1 < ./rep1_lib.txt
#readarray LIB_REP2 < ./rep2_lib.txt
#readarray LIB_CONTROL1 < ./control1_lib.txt
#readarray LIB_CONTROL2 < ./control2_lib.txt
#
#LIB_REP1=( $(echo ${LIB_REP1[@]} | tr " " "\n" | tr "\n" " "))
#LIB_REP2=( $(echo ${LIB_REP2[@]} | tr " " "\n" | tr "\n" " "))
#LIB_CONTROL1=( $(echo ${LIB_CONTROL1[@]} | tr " " "\n" | tr "\n" " "))
#LIB_CONTROL2=( $(echo ${LIB_CONTROL2[@]} | tr " " "\n" | tr "\n" " "))


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


for i in "${!LIB_REP1[@]}"; do

    R1=${LIB_REP1[$i]}
    R2=${LIB_REP2[$i]}
	IDR_DIR=$BASE_DIR/IDR_${R1}_${R2}
	echo -e "\nidr.${R1}_${R2}.truerep.out : name of log file"
	echo -e "\ndir for logfile : $IDR_DIR/idr.${R1}_${R2}.truerep.out" 
	echo -e "\ncheck all idr peaks agnostic of cutoff: $IDR_DIR/${R1}.filt.nodup.srt.SE_VS_${R2}.filt.nodup.srt.SE.IDR0.narrowPeak.gz\n"
	echo -e "\ncheck idr passed peaks with cutoff 0.05 (Real peaks for the analysis) : $IDR_DIR/${R1}.filt.nodup.srt.SE_VS_${R2}.filt.nodup.srt.SE.IDR0.02.filt.narrowPeak.gz\n"

	# True Rep IDR
	bsub $BSUB_OPTIONS -o "$IDR_DIR/idr05.${R1}_${R2}.truerep.out" -J "IDR $R1 $R2 true reps" ./Idr.main_05.sh $IDR_DIR/$R1.filt.nodup.srt.SE.narrowPeak.gz $IDR_DIR/$R2.filt.nodup.srt.SE.narrowPeak.gz $IDR_DIR/${R1}_${R2}.Rep0.narrowPeak.gz signal.value $IDR_DIR

	# Rep 1 self-pseudoreps (Rep1 psuedorep1 vs Rep1 pseudorep2)
	bsub $BSUB_OPTIONS -o "$IDR_DIR/idr05.${R1}_${R2}.rep1_pseudo.out" -J "IDR $R1 $R2 rep1 pseudo" ./Idr.main_05.sh $IDR_DIR/$R1.filt.nodup.SE.pr1.narrowPeak.gz $IDR_DIR/$R1.filt.nodup.SE.pr2.narrowPeak.gz $IDR_DIR/$R1.filt.nodup.srt.SE.narrowPeak.gz signal.value $IDR_DIR

	# Rep 2 self-pseudoreps (Rep2 pseudorep1 vs Rep2 pseudorep2) 
	bsub $BSUB_OPTIONS -o "$IDR_DIR/idr05.${R1}_${R2}.rep2_pseudo.out" -J "IDR $R1 $R2 rep2 pseudo" ./Idr.main_05.sh $IDR_DIR/$R2.filt.nodup.SE.pr1.narrowPeak.gz $IDR_DIR/$R2.filt.nodup.SE.pr2.narrowPeak.gz $IDR_DIR/$R2.filt.nodup.srt.SE.narrowPeak.gz signal.value $IDR_DIR

	# Pooled psuedoreps ( Pooled pseudorep1 (PR1 of Rep1 + PR1 of Rep2) vs Pooled true tag align file(Rep1_TA_file + Rep2_TA_file)
	bsub $BSUB_OPTIONS -o "$IDR_DIR/idr05.${R1}_${R2}.pool_pseudo.out" -J "IDR $R1 $R2 pooled pseudo" ./Idr.main_05.sh $IDR_DIR/${R1}_${R2}.Rep0.pr1.narrowPeak.gz $IDR_DIR/${R1}_${R2}.Rep0.pr2.narrowPeak.gz $IDR_DIR/${R1}_${R2}.Rep0.narrowPeak.gz signal.value $IDR_DIR

done



# Pooled psuedoreps ( Pooled pseudorep2 (PR2 of Rep1 + PR2 of Rep2) vs Pooled true tag align file(Rep1_TA_file + Rep2_TA_file)
# Pooled true Tag align (Rep1_TA_file + Rep2_TA_file)  vs Pooled control tag align (Control1_TA_file _ Control2_TA_file)
