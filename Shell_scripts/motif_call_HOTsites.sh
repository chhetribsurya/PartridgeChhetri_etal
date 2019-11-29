#!/bin/bash

file_path=$1
MEME_DB_FILE=$2
meme_target_dir=$3
TRUE_HOT_FILE=$4

if [[ ! -d $meme_target_dir ]]; then mkdir -p $meme_target_dir;fi

export BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin"
export MEME_SUITE_PATH="/gpfs/gpfs2/software/meme-4.11.4/bin"
#export GENOME="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/hg38_refgenome/hg38-male.fa"
export GENOME="/gpfs/gpfs1/home/schhetri/encode_k562_data_analysis/bin/ref_fasta/male.hg19.fa"
#export TRUE_HOT_FILE="MasterBedFile.${sample}.bed.srt.mrg.filt.70.hot.final"

cat "$file_path" | while read each;do
        #echo "debug $each"
        tf_name=$(echo $each | cut -f1 -d " ")
        #echo "check tf_name : $tf_name"
        cisbp_id=$(echo $each | cut -f2 -d " ")
        $MEME_SUITE_PATH/meme-get-motif -id $cisbp_id $MEME_DB_FILE > ${meme_target_dir}/${tf_name}_${cisbp_id}.meme
done

combined_meme_file=$(basename $file_path .txt).combined.meme
cat $(ls $meme_target_dir/*.meme | grep -v "combined.meme") > $meme_target_dir/$combined_meme_file
#cat $meme_target_dir/*.meme > $meme_target_dir/$combined_meme_file

# Conversion of HOT regions to fasta:
hot_file=$(dirname $file_path)/$TRUE_HOT_FILE
hot_file=$(ls $hot_file)
echo "Hot file :$hot_file"

export ORIGINAL_HOT_FASTA_FILE=${hot_file}.fasta
export FIMO_DIR=$(dirname $file_path)/fimo_dir

if [[ ! -d $FIMO_DIR ]];then mkdir $FIMO_DIR;fi

# Fasta extraction of HOT region:
$BEDTOOLS_PATH/fastaFromBed -fi $GENOME -bed $hot_file -fo $ORIGINAL_HOT_FASTA_FILE
echo -e "\nFasta extraction completed"

# Find motif occurrence of given permuted tfs at HOT regions:
echo "Hot fasta file : $ORIGINAL_HOT_FASTA_FILE"
$MEME_SUITE_PATH/fimo --parse-genomic-coord --verbosity 1 --oc $FIMO_DIR --thresh 1e-6 $meme_target_dir/$combined_meme_file ${ORIGINAL_HOT_FASTA_FILE}

echo -e "\nFimo dir for HOT site: $FIMO_DIR"
echo -e "\nFimo task completed!"

index_id=$(basename $hot_file | awk -F. '{print $(NF-1)}')
hotmotif_file=$(dirname $file_path)/Hotsites.${index_id}.fimo.txt

# Merging of motif to avoid adjacent motif region being overrepresented:
awk '{print $3,$4,$5,$8,$9,$2}' OFS="\t" $FIMO_DIR/fimo.txt | tail -n+2 > $hotmotif_file
sort -k1,1 -k2,2n $hotmotif_file > ${hotmotif_file}.srt
$BEDTOOLS_PATH/mergeBed -i ${hotmotif_file}.srt -c 6,6,6 -o count,collapse,count_distinct > ${hotmotif_file}.srt.mrg

# Filter hot regions to keep analysis cleaner:
awk '{print $1,$2,$3,$5,$6}' OFS="\t" ${hot_file} > ${hot_file}.final

$BEDTOOLS_PATH/intersectBed -a ${hot_file}.final -b ${hotmotif_file}.srt.mrg -wa -wb > $(dirname $file_path)/Hot_fimomotif_intersect.${index_id}.bed.txt
$BEDTOOLS_PATH/intersectBed -a ${hot_file}.final -b ${hotmotif_file}.srt.mrg -v > $(dirname $file_path)/Hot_fimomotif_outersect.${index_id}.bed.txt
echo -e "\nIntersection of Hot sites and motif call completed!"
echo -e "Final file : $(dirname $file_path)/Hot_fimomotif_intersect.${index_id}.bed.txt"





