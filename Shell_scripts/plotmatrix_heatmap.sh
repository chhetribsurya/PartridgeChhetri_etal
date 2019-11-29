# run compute matrix to collect the data needed for plotting

bigwig_dir="/gpfs/gpfs1/home/schhetri/paper_revision_3/deeptools/nurd_complex"
output_dir="/gpfs/gpfs1/home/schhetri/paper_revision_3/deeptools/nurd_complex/deeptools_matrix_plots_final"
output_file="matrix_Gatad2a_signal_FOXA3only.gz"
output_tab_file=$(basename "matrix_Gatad2a_signal_FOXA3only.gz" .gz).tab
output_plot_file=$(basename "matrix_Gatad2a_signal_FOXA3only.gz" .gz).png

if [[ ! -d $output_dir ]];then mkdir -p $output_dir; fi

computeMatrix reference-point --referencePoint center -p max -b 300 -a 300 \
                                                        \
                                -R $bigwig_dir/FOXA3_Only.bed \
                                -R $bigwig_dir/GATAD2A_Only.bed \
                                                            \
                                -S  $bigwig_dir/GATAD2A_FLAG.bw \
                                    $bigwig_dir/FOXA3_FLAG.bw  \
                                    $bigwig_dir/SMAD4_iso1_FLAG.bw \
                                    $bigwig_dir/RARA_FLAG.bw   \
                                    $bigwig_dir/SOX13_iso1_FLAG.bw \
                                    $bigwig_dir/ZNF219_FLAG.bw \
                                    $bigwig_dir/ARID5B_FLAG.bw \
                                            \
                                --skipZeros \
                                --numberOfProcessors max \
                                -o $output_dir/$output_file \
                                --outFileNameMatrix $output_dir/$output_tab_file
                              



 plotHeatmap -m $output_dir/$output_file \
     -o $output_dir/$output_plot_file \
     --regionsLabel "FOXA3 only GATAD2A only" \
     --refPointLabel "Center" \
     --xAxisLabel "Peak Distance(bp)" \
     #--whatToShow 'heatmap and colorbar' \
     #--zMin -3 --zMax 3 
     #--colorMap RdBu \
     #--kmeans 4 \
     #--boxAroundHeatmaps no \
