#!/bin/bash

sample="HepG2New"
file_path=$1
MEME_DB_FILE=$2
meme_target_dir=$3

export BEDTOOLS_PATH="/gpfs/gpfs2/software/bedtools2-2.20.0/bin"
export MEME_SUITE_PATH="/gpfs/gpfs2/software/meme-4.11.4/bin"
#export GENOME="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/hg38_refgenome/hg38-male.fa"
export GENOME="/gpfs/gpfs1/home/schhetri/encode_k562_data_analysis/bin/ref_fasta/male.hg19.fa"
#export TRUE_HOT_FILE="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/cisbpmotif_id_extract/MasterBedFile.HepG2.bed.srt.mrg.filt.78.hot.final"
export TRUE_HOT_FILE="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill/$sample/cisbpmotif_id_extract/MasterBedFile.HepG2New.bed.srt.mrg.filt.70.hot.final"

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
hot_file=$(dirname $file_path)/MasterBedFile*.hot
index_id_temp=$(basename $hot_file | cut -f3-4,5-6 -d ".")
hot_file=$(ls ${hot_file})
$BEDTOOLS_PATH/intersectBed -a ${hot_file} -b ${TRUE_HOT_FILE} -u > ${hot_file}sub

$BEDTOOLS_PATH/intersectBed -a ${TRUE_HOT_FILE} -b ${hot_file} -u > ${hot_file}.fracOverlap
total_truehot_count=$(wc -l < ${TRUE_HOT_FILE} | sed 's/^ *//g')
total_hot_overlap=$(wc -l < ${hot_file}.fracOverlap | sed 's/^ *//g')
percent_overlap=$(echo "scale=4;${total_hot_overlap}/${total_truehot_count}*100" | bc)
echo -e "${index_id_temp}\nTotal_hot_count/Total_truehot_count : $total_hot_overlap/$total_hot_overlap*100 = $percent_overlap%"
echo -e "${index_id_temp}\t${total_hot_overlap}/${total_truehot_count}*100\t$percent_overlap" >> $(dirname ${TRUE_HOT_FILE})/TRUEHOT_SUBHOT.fracOverlap.txt

hot_file=$(ls ${hot_file}sub)
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

index_id=$(basename $hot_file | cut -f3-4,5-6 -d ".")
hotmotif_file=$(dirname $file_path)/Hotsites.${index_id}.fimo.txt

# Merging of motif to avoid adjacent motif region being overrepresented:
awk '{print $3,$4,$5,$8,$9,$2}' OFS="\t" $FIMO_DIR/fimo.txt | tail -n+2 > $hotmotif_file
sort -k1,1 -k2,2n $hotmotif_file > ${hotmotif_file}.srt
$BEDTOOLS_PATH/mergeBed -i ${hotmotif_file}.srt -c 6,6,6 -o count,collapse,count_distinct > ${hotmotif_file}.srt.mrg

# Filter hot regions to keep analysis cleaner:
awk '{print $1,$2,$3,$6}' OFS="\t" ${hot_file} > ${hot_file}.final

$BEDTOOLS_PATH/intersectBed -a ${hot_file}.final -b ${hotmotif_file}.srt.mrg -wa -wb > $(dirname $file_path)/Hot_fimomotif_intersect.${index_id}.bed.redo.txt
echo -e "\nIntersection of Hot sites and motif call completed!"
echo -e "Final file : $(dirname $file_path)/Hot_fimomotif_intersect.${index_id}.bed.redo.txt"





