#!/bin/bash

$MEME_PATH/tomtom -oc ${TOMTOM_OUTPUT_DIR} -png -eps -min-overlap 5 -dist pearson -evalue -thresh 1 -no-ssc -bfile ${BG_FILE} ${MEME_FILE} ${ALL_TARGET_FILES}
