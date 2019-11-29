#!/bin/bash

RUN_PATH=$(pwd)
SAMPLE="HepG2full"
#SAMPLE="HepG2"
#SAMPLE="K562"
#SAMPLE="GM12878"
#SAMPLE="H1hESC"

UNIQ_TF_DIR="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/$SAMPLE/bedfiles"  #ENCODE_K562.IDR0.02.filt.narrowPeak_MTA3.fasta
OUTPUT_DIR="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/$SAMPLE"
BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin"
#REFERENCE_GENOME_PATH="/gpfs/gpfs1/home/schhetri/encode_k562_data_analysis/bin/ref_fasta"

#MERGED_TFs_DIR=${UNIQ_TF_DIR}/merged_TFs_dir
if [[ ! -d ${OUTPUT_DIR} ]];then
    mkdir -p ${OUTPUT_DIR}
fi

MERGED_TFs_DIR=${OUTPUT_DIR}/merged_TFs_dir_final
if [[ ! -d ${MERGED_TFs_DIR} ]];then
    mkdir -p ${MERGED_TFs_DIR}
fi
MASTER_BED_FILE=${MERGED_TFs_DIR}/MasterBedFile.${SAMPLE}.bed
MASTER_BED_FILE_SRT=${MASTER_BED_FILE}.srt
MASTER_BED_FILE_SRT_MRG=${MASTER_BED_FILE_SRT}.mrg
MASTER_BED_FILE_SRT_MRG_FA=${MASTER_BED_FILE_SRT_MRG}.fa
MASTER_BED_FILE_SRT_MRG_FILT=${MASTER_BED_FILE_SRT_MRG}.filt
MASTER_BED_FILE_SRT_MRG_FILT_FA=${MASTER_BED_FILE_SRT_MRG}.filt.fa

printf "\nProcessing and Combining TF files from %s ...\n" ${UNIQ_TF_DIR}
### Printing 2 extra columns with {TFname}_{linenum} and {TFname} - useful for retaining merged peaks identity:
for each in $(ls ${UNIQ_TF_DIR}/*narrowPeak*); do
    #TF_name=$( echo $each | awk -F '[narrowPeak]' '{print $NF}' | cut -f2- -d "_" ); 
    TF_name=$(echo $each|sed -rn 's/.*narrowPeak_(.*)/\1/p'| sed 's/\[FLAG\]//'); #echo "Cool : $TF_name"
    awk -v TF="$TF_name" '{print $0"\t"TF "_" NR"\t"TF}' $each ; 
done > ${MASTER_BED_FILE}
#
printf "Now sorting masterbed file...\n"
### Sort all bed files for bedtools-mergeBed to work on:
sort -k1,1 -k2,2n ${MASTER_BED_FILE} > ${MASTER_BED_FILE_SRT}

### Merge all the sorted peaks of TF binding sites for merged hot and non-hot sites generation:
printf "Now Merging sorted masterbed file...\n"
${BEDTOOLS_PATH}/mergeBed -i ${MASTER_BED_FILE_SRT} -c 11,11,12 -o count,collapse,count_distinct > ${MASTER_BED_FILE_SRT_MRG}

### Filter all merged peaks spreading more than 2kb regions:
awk 'BEGIN{OFS="\t"} {chromStart=$2; chromEnd=$3; width=chromEnd-chromStart; if(width<=2000) print $0}' ${MASTER_BED_FILE_SRT_MRG} > ${MASTER_BED_FILE_SRT_MRG_FILT}

### Generate fasta for merged TF files for SVM prediction using model weights using SVM:
#printf "Now generating fasta for merged masterbed file...\n"
#${BEDTOOLS_PATH}/fastaFromBed -fo ${MASTER_BED_FILE_SRT_MRG_FA} -fi ${REFERENCE_GENOME_PATH}/male.hg19.fa -bed ${MASTER_BED_FILE_SRT_MRG}
#${BEDTOOLS_PATH}/fastaFromBed -fo ${MASTER_BED_FILE_SRT_MRG_FILT_FA} -fi ${REFERENCE_GENOME_PATH}/male.hg19.fa -bed ${MASTER_BED_FILE_SRT_MRG_FILT}
echo -e "Merged MasterBedFile Generated!!\n"
