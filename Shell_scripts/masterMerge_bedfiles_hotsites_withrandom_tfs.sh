#!/bin/bash

RUN_PATH=$(pwd)
tf_count=$1
random_seed=$2
new_output_dir=$3
SAMPLE=$4

UNIQ_TF_DIR="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/$SAMPLE/bedfiles"  #ENCODE_K562.IDR0.02.filt.narrowPeak_MTA3.fasta
OUTPUT_DIR=$new_output_dir
BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin"
#REFERENCE_GENOME_PATH="/gpfs/gpfs1/home/schhetri/encode_k562_data_analysis/bin/ref_fasta"

MASTER_BED_FILE=${OUTPUT_DIR}/MasterBedFile.${SAMPLE}.tfcount.${tf_count}.perm.${random_seed}.bed
MASTER_BED_FILE_SRT=${MASTER_BED_FILE}.srt
MASTER_BED_FILE_SRT_MRG=${MASTER_BED_FILE_SRT}.mrg
MASTER_BED_FILE_SRT_MRG_FA=${MASTER_BED_FILE_SRT_MRG}.fa
MASTER_BED_FILE_SRT_MRG_FILT=${MASTER_BED_FILE_SRT_MRG}.filt
HOT_BED_FILE_SRT_MRG_FILT=${MASTER_BED_FILE_SRT_MRG}.filt.hot
MASTER_BED_FILE_SRT_MRG_FILT_FA=${MASTER_BED_FILE_SRT_MRG}.filt.fa

printf "\nProcessing and Combining TF files for Random tf from %s ...\n" ${UNIQ_TF_DIR}

get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}

### Printing 2 extra columns with {TFname}_{linenum} and {TFname} - useful for retaining merged peaks identity:
randomcount_tf_list=$(ls ${UNIQ_TF_DIR}/*narrowPeak* | shuf --random-source=<(get_seeded_random ${random_seed}) -n ${tf_count})
echo $randomcount_tf_list | tr " " "\n" > $OUTPUT_DIR/tfcount.${tf_count}.perm.${random_seed}.filelist.txt

for each in ${randomcount_tf_list}; do
    TF_name=$(echo $each|sed -rn 's/.*narrowPeak_(.*)/\1/p'| sed 's/\[FLAG\]//'); #echo "Cool : $TF_name"
    awk -v TF="$TF_name" '{print $0"\t"TF "_" NR"\t"TF}' $each ; 
done > ${MASTER_BED_FILE}

printf "Now sorting masterbed file...\n"
### Sort all bed files for bedtools-mergeBed to work on:
sort -k1,1 -k2,2n ${MASTER_BED_FILE} > ${MASTER_BED_FILE_SRT}

### Merge all the sorted peaks of TF binding sites for merged hot and non-hot sites generation:
printf "Now Merging sorted masterbed file...\n"
${BEDTOOLS_PATH}/mergeBed -i ${MASTER_BED_FILE_SRT} -c 11,11,12 -o count,collapse,count_distinct > ${MASTER_BED_FILE_SRT_MRG}

### Filter all merged peaks spreading more than 2kb regions:
awk 'BEGIN{OFS="\t"} {chromStart=$2; chromEnd=$3; width=chromEnd-chromStart; if(width<=2000) print $0}' ${MASTER_BED_FILE_SRT_MRG} > ${MASTER_BED_FILE_SRT_MRG_FILT}
echo -e "Merged MasterBedFile Generated!!\n"

hotsite_cutoff=$(echo $tf_count | awk '{print $1/2}')
echo -e "\nHOT site cutoff set = $hotsite_cutoff TFs"
awk -v var=$hotsite_cutoff '{if($6 >= var) print $0}' $MASTER_BED_FILE_SRT_MRG_FILT > $HOT_BED_FILE_SRT_MRG_FILT

### Generate fasta for merged TF files for SVM prediction using model weights using SVM:
#printf "Now generating fasta for merged masterbed file...\n"
#${BEDTOOLS_PATH}/fastaFromBed -fo ${MASTER_BED_FILE_SRT_MRG_FA} -fi ${REFERENCE_GENOME_PATH}/male.hg19.fa -bed ${MASTER_BED_FILE_SRT_MRG}
#${BEDTOOLS_PATH}/fastaFromBed -fo ${MASTER_BED_FILE_SRT_MRG_FILT_FA} -fi ${REFERENCE_GENOME_PATH}/male.hg19.fa -bed ${MASTER_BED_FILE_SRT_MRG_FILT}

