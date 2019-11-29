#!/usr/bin/bash

#using hg19 ref seq
sample="HepG2full"
main_dir="/gpfs/gpfs1/home/schhetri/for_ENCODE_jill"
merged_tf_file="$main_dir/$sample/merged_TFs_dir_final/MasterBedFile.${sample}.bed.srt.mrg.filt"

### Set the path for all the tools to be consistent:
export MEME_DB_FILE="/gpfs/gpfs1/home/schhetri/Tools/meme_4.11.2/db/motif_databases/motif_databases/CIS-BP/Homo_sapiens_edit.meme"
#export MEME_DB_FILE="/gpfs/gpfs1/home/schhetri/for_jacob/quartile_motif_analysis/Homosapiens_combined_database_CISBP_JASPAR_edit.meme"
#export TRUE_HOT_FILE="MasterBedFile.${sample}.bed.srt.mrg.filt.70.hot.final"

### Filter hot file >=X/2 factors; here >=104 factors;
awk '{if($6>=70) print $0}' $merged_tf_file > "$main_dir/$sample/cisbpmotif_id_extract/$(basename $merged_tf_file).hot"

export TRUE_HOT_FILE="$(basename $merged_tf_file).hot"


### Set the pat for cisbpID containing TFs:
tf_cismemeid_path="$main_dir/$sample/cisbpmotif_id_extract/tf_list_withCisBP_id.txt"

target_dir=$(dirname $tf_cismemeid_path)
file_name=$(basename $tf_cismemeid_path)
meme_target_dir=$target_dir/all_nonpermuted_tf_meme_pwms
if [[ ! -d $meme_target_dir ]];then mkdir -p $meme_target_dir;fi 
bsub -We 2:00 -q c7priority -J "motif extraction for $file_name" -o ./HOTsite_only_getmotif_step.log ./motif_call_HOTsites.sh $tf_cismemeid_path $MEME_DB_FILE $meme_target_dir ${TRUE_HOT_FILE}
#    done





 
