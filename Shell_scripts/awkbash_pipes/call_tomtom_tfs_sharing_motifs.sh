#!/bin/bash

RUN_PATH=`pwd`
export MEME_FILE_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_chip_motif_analysis_total/unique_TFs" 
export OUTPUT_DIR="/gpfs/gpfs1/home/schhetri/for_chris/batch_I/meme_motif_sharing_bw_tfs"
export MEME_PATH="/gpfs/gpfs2/software/meme-4.11.4/bin/"

for each_file in $(ls -d $MEME_FILE_DIR/*narrowPeak*); do 
    TF_NAME=$(basename $each_file | sed -rn 's/.*_narrowPeak_(.*$)/\1/p')
    TOMTOM_OUTPUT_DIR=$OUTPUT_DIR/$TF_NAME
    echo -e "\nProcessing: $TF_NAME"

    if [[ ! -d $TOMTOM_OUTPUT_DIR ]];then
        mkdir -p $TOMTOM_OUTPUT_DIR
    fi

    ORIGINAL_MOTIF_FILE=$each_file/meme_chip/meme_out/meme.xml
    MEME_FILE=$TOMTOM_OUTPUT_DIR/pssm_file_${TF_NAME}.meme

    $MEME_PATH/meme-get-motif -all $ORIGINAL_MOTIF_FILE | sed -r "s/^MOTIF.*([[:digit:]])/MOTIF ${TF_NAME}.\1 MEME/" > $MEME_FILE 
    echo "Task completed for ${TF_NAME}!!" 
done

### Once all the pssm meme file is formed; convert all those meme files to target or database files:
echo -e "\n\nMerging all the meme files extracted in the last for loop stage...\n"
cat $(ls -d $OUTPUT_DIR/*/pssm_file_*.meme | xargs) > $OUTPUT_DIR/custom_TF_motif_concat_db.meme

### Remove all the repetitions of mem_version, alphabet lines and so on..."
$RUN_PATH/remove_multiple_alphabets.awk $OUTPUT_DIR/custom_TF_motif_concat_db.meme > $OUTPUT_DIR/custom_TF_motif_edited_db.meme

### Replace background alphabet frequency with uniform ACGT frequency i.e with 0.25 for each"
sed -r 's/ A [0-9]\.[[:digit:]]+ C [0-9]\.[[:digit:]]+ G [0-9]\.[[:digit:]]+ T [0-9]\.[[:digit:]]+/ A 0.25000 C 0.25000 G 0.25000 T 0.25000/' $OUTPUT_DIR/custom_TF_motif_edited_db.meme > $OUTPUT_DIR/custom_TF_motif_final_db.meme
echo -e "\nCustom motif database generated...."

### Copy all the background model file to respective TF dir to sort brackets issue in the filenames(case: BG_FILE=$each_file/meme_chip/*.model)
for each_model in $(ls $MEME_FILE_DIR/*/meme_chip/*.model);do
    TF_name=$(echo $each | sed -rn 's/.*_narrowPeak_(.*)_background_.*/\1/p')
    cp $each_model $TOMTOM_OUTPUT_DIR/${TF_name}/${TF_name}_background.model 
done

### Running Tomtom algorithm for quantifying the similarity between the motifs for each TFs:
export ALL_TARGET_FILES=$(ls -d $OUTPUT_DIR/*/pssm_file_*.meme | xargs)
#BSUB_OPTIONS="-We 24:00 -q c7normal [mem=8096]" # Using new cluster
export BSUB_OPTIONS="-We 24:00 -q c7normal -R rusage[mem=10000]" # Using new cluster

for each_file in $(ls -d $MEME_FILE_DIR/*narrowPeak*); do 
    export TF_NAME=$(basename $each_file | sed -rn 's/.*_narrowPeak_(.*$)/\1/p')
    export TOMTOM_OUTPUT_DIR=$OUTPUT_DIR/$TF_NAME
    echo -e "\n\nTomtom running against all target files for: $TF_NAME\n"

    export ORIGINAL_MOTIF_FILE=$each_file/meme_chip/meme_out/meme.xml
    export BG_FILE=$TOMTOM_OUTPUT_DIR/${TF_NAME}_background.model
    export MEME_FILE=$TOMTOM_OUTPUT_DIR/pssm_file_${TF_NAME}.meme

    bsub $BSUB_OPTIONS -J "Motif comparision with other TFs motif db" -o ${OUTPUT_DIR}/tomtom_run_motif_comparison.out $RUN_PATH/tomtom_tfs_sharing_motifs_call.sh
    #$MEME_PATH/tomtom -oc ${TOMTOM_OUTPUT_DIR} -png -eps -min-overlap 5 -dist pearson -evalue -thresh 1 -no-ssc -bfile ${BG_FILE} ${MEME_FILE} ${ALL_TARGET_FILES}
    echo -e "\nTask completed for ${TF_NAME}!!" 
done

echo -e "\n\n\nJobs successfully completed!!!\n"
