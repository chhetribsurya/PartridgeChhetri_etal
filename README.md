# Occupancy maps of 208 chromatin-associated proteins in one human cell type

## All the codes used in this project, and study, are split into their language specific platforms : Python, R and Shell (bash, awk)

##### `Code execution may require installation of following softwares and dependencies:`

```
   - Python 2.7
   - Python >= 3.5
   - R >= 3.5
   - Shell Env (unix, bash) 
   - boost C++ libraries
   - bwa >= 0.7.12
   - bowtie2 >= 2.1.0
   - spp >= 1.10.1
   - idr >= 2.0.2
   - meme-4.11.4
   - centrimo >= 4.11.4
   - tomtom >= 4.11.4
   - samtools >= 1.3
   - bedtools2 >= 2.20.0
   - fimo >= 4.11.4
   - phantompeakqualtools >= 2.0
   - picard-tools >= 1.88 
   - trim_galore >=0.3.7 
   - cutadapt >= 1.16
   - fastqc >= 0.10.1
   - deeptools >= 3.0.0
   - kmersvm; gkmSVM >= 2.0
   - pybedtools >= 0.7.10
   - pandas >= 0.20.3
   - numpy >= 1.14.0
   - scipy >= 0.19.1
   - scikit-learn >= 0.19.0 
   - ggplot2 >= 3.1.0; gplots >= 3.0.1
   - dplyr >= 0.7.5; gtools >= 3.0.1
   - ranger >= 0.10.1
   - ComplexHeatmap >= 1.18.0
   - circlize >= 0.4.3
   - igraph >= 1.2.1
```

##### `Databases used in the study:`

```
   cisbp
   jaspar
   meme
```

##### `Inputs for the study are chiefly TFBS (Transcription factor binding sites) as a peak bed files mapped by ChIP-Seq and CETCh-Seq assay`


```
   Some analysis implements the standalone scripts, and some the combination of multiple or chaining scripts, including pipes.

   Codes built here were not originally intended for packages, but serves as a framework, and maybe considered 

   in future for more extensive analysis purpose, as required.
```

`Note : Please cite the orginal softwares and programs if used in the study, for example, SPP, IDR, Bedtools, meme-suite, phantompeakqualtools and kmer/gkm-SVM... based papers, to which some of the codes here are adapted from as well. The relevant references are here below:`


1. Landt SG, Marinov GK, Kundaje A, et al. ChIP-seq guidelines and practices of the ENCODE and modENCODE consortia. Genome Res. 2012;22(9):1813-1831. doi:10.1101/gr.136184.111

2. Li QH, Brown JB, Huang HY, Bickel PJ. Measuring reproducibility of high-throughput experiments. Ann. Appl. Stat. 2011; 5(3):1752-1779.

3. Kharchenko PV, Tolstorukov MY, Park PJ. Design and analysis of ChIP-seq experiments for DNA-binding proteins. Nat Biotechnol. 2008 Dec;26(12):1351-9. PMID: 19029915; PMCID: PMC2597701

4. Dale, R. K., Pedersen, B. S. & Quinlan, A. R. Pybedtools: a flexible Python library for manipulating genomic datasets and annotations. Bioinformatics 27, 3423–3424 (2011).

5. Quinlan, A. R. & Hall, I. M. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics 26, 841–842 (2010).

6. Fletez-Brant, C., Lee, D., McCallion, A. S. & Beer, M. A. kmer-SVM: a web server for identifying predictive regulatory sequence features in genomic data sets. Nucleic Acids Res. 41, W544–W556 (2013).

7. Machanick, P. & Bailey, T. L. MEME-ChIP: motif analysis of large DNA datasets. Bioinformatics 27, 1696–1697 (2011).

8. Ma, W., Noble, W. S. & Bailey, T. L. Motif-based analysis of large nucleotide data sets using MEME-ChIP. Nat. Protocols 9, 1428–1450 (2014).

9. Ghandi, M.†, Lee, D.†, Mohammad-Noori, M. & Beer, M. A. Enhanced Regulatory Sequence Prediction Using Gapped k-mer Features. PLoS Comput Biol 10, e1003711 (2014). doi:10.1371/journal.pcbi.1003711 † Co-first authors

10. Lee, D. LS-GKM: A new gkm-SVM for large-scale Datasets. Bioinformatics btw142 (2016). doi:10.1093/bioinformatics/btw142

11. Partridge EC, Chhetri SB, Prokop JW, et al. Occupancy maps of 208 chromatin-associated proteins in one human cell type. Nature. 2020;583(7818):720-728. doi:10.1038/s41586-020-2023-4

