#!/usr/bin/env Rscript

# Install packages if not installed prior:
list.of.packages <- c("data.table", "dplyr", "ComplexHeatmap", "circlize", "colorspace", "GetoptLong", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(colorspace)
library(GetoptLong)
library(optparse)


args <- commandArgs(trailingOnly=TRUE)

# Check if there is at least one argument: if not return an error
option_list <- list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
    make_option(c("-o", "--out"), type="character", default="Ideas_state_piechart_heatmap_figure.pdf", 
              help="output file name [default= %default]", metavar="character")
    make_option(c("-outdir", "--output dir"), type="character", default=".", 
              help="output dir name [default= %default]", metavar="character")
); 
 
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);


if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Load TF annoation file:
# anno_df <- fread("~/Dropbox/misc/TFs_Annotation_file.txt", sep="\t")
# anno_df <- anno_df %>% mutate(label2=ifelse(Category=="DBF", "black", "grey"))

# Load concatenated file (previous output):
input_file <- args[1]

# Read input files:
read_file <- fread(input_file, sep="\t", header=TRUE)
# read_file_merged <- merge(read_file, anno_df, by.x="tf_name", by.y="Target", all.x=TRUE)
# read_file_merged$tf_name_edit <- paste(read_file_merged$tf_name, read_file_merged$annotation)

read_df <- as.data.frame(read_file)
data <- read_df[,2:ncol(read_df)] # Subtract last 4 columns
mat_data <- as.matrix(data)

# Binding variance across the cis-regulatory region for any given factor
mat_data <- t(apply(mat_data, 1, scale))
colnames(mat_data) <- colnames(data)

rownames(mat_data) <- read_df$tf_name
rownames(mat_data)
colnames(mat_data)

# Input output dir target:
output_dir <- "." # current dir
output_file_name <- file.path(args[3], args[2])       
pdf(output_file_name)

# Build heatmap
Heatmap(mat_data, name="Binding Z-score", 
    column_title="IDEAS Genomic States",
    row_title="Transcription Factors",
    row_title_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_gp = gpar(fontsize = 1.4),
    column_names_gp = gpar(fontsize = 6)
) 

dev.off()
