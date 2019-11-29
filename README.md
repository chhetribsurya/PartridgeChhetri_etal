# "Occupancy maps of 208 chromatin-associated proteins in one human cell type"

## All the codes used in this project, and study, are segregated into their language specific platforms : Python, R and Shell (bash,awk)

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

`Note : Please cite the orginal software and program if used in the study, for example, SPP, IDR, Bedtools, meme-suite and kmer/gkm-SVM... based papers`

