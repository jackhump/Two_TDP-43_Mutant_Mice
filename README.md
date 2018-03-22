# Fratta, Sivakumar, Humphrey, Lo, Ricketts, et al, 2018
# **Endogenous TDP-43 mutant mice have novel gain of splicing function and ALS characteristics in vivo**
# in press at EMBO Journal
# 2018

This GitHub repository contains all scripts to recreate the key bioinformatic analyses of the paper.

all scripts written by Jack Humphrey, Kitty Lo, Prasanth Sivakumar and Shannon Edwards

## Dependencies:

Mac OS X or Linux

bedtools v2.24.0
bigWigSummary (http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigSummary)

*R >= 3.3.2*
- dplyr
- ggplot2
- stringr
- biostrings
- data.table
- gridExtra
- optparse
- smoother (no longer available on CRAN for R 3.4.3, have to build from tarball - sorry! )
- TTR

*python 2.7.10*
- pybedtools
- argparse


## Installation

1. clone github repository to your computer
1. unzip the data:

```
tar -xzvf data.tar.gz
```

3. download the PhyloP 60way conservation scores from UCSC:

```
wget http://hgdownload.cse.ucsc.edu/goldenpath/mm10/phyloP60way/mm10.60way.phyloP60way.bw
```


## A note on data pre-processing

All RNA-seq data was processed using the Plagnol lab RNA-seq pipeline (github.com/plagnollab/RNASeq_pipeline). 

Samples were aligned to mm10 with STAR.

Reads that overlap exons were counted with HTSeq. 

Differential gene expression was computed with DESeq2.

Differential splicing was computed with SGSeq.

All iCLIP data was downloaded from iCOUNT in low FDR clusters.

## Analysis scripts used in the paper:

###  - 5A,B) RNA maps 	

-  `RNA_maps/whole_intron_cluster_coverage.py`	
-  `RNA_maps/plot_whole_intron_coverage.R`
-  `RNA_maps/create_RNAmaps.sh`

`whole_intron_cluster_coverage.py` takes a BED file of exons or introns and intersects them with a second BED file of iCLIP clusters.

`plot_whole_intron_coverage` composes the RNAmaps from a set of introns and exons. Per-gene coverages are normalised by the total number of genes to create a per-nucleotide coverage.

Reproduce the plots in 5A) and 5B) by running: 

```
cd RNAmaps/
sh create_RNAmaps.sh
```


### 5C) conservation
- `conservation/conservation.R`

Calculates per-exon conservation score for groups of exons.

bigWigSummary must be downloaded from http://hgdownload.soe.ucsc.edu/admin/exe/macOSX.x86_64/bigWigSummary and installed on your $PATH

Reproduce figure 5C):

```
cd conservation/
Rscript conservation.R
```

### 5D) Differential expression 
- `differential_expression/skiptic_cryptic_expression.R`

Running this script compares the DESeq2 differential gene expression results with the lists of genes that contain extreme splicing events.

### 5E) functional prediction
- `NMD_prediction/NMD_prediction.R`
- `NMD_prediction/bring_everything_together.R`

This needs explaining. 

### Expanded View 1
A,D) pie charts
- `differential_splicing/makePieChartsAllEvents.R`

creates pie charts of splicing variant types. 


### Expanded View 4
A,B) long genes
- `differential_expression/all_mice_TDP-43_long_genes.R`

Takes each DESeq2 differential expression results file and creates a signed Z-score for each gene ( Z = qnorm(1 - (P / 2)  ) with the sign of the direction of the fold change. Genes are then ranked and binned into groups of 200 genes. The mean intron length for each gene according to GENCODE annotation is then matched and compared.

### Appendix S6 - MA plots 
- `ma_plots/make_MA_plots.R`

Takes the SGSeq differential splicing results file and plots the mean reads that support each splice event against the log2 fold change of that splice event.


### Appendix S7 - Permutation testing
- permutation/permute_qq.R
- permutation/permute_100_sgseq.sh
- permutation/plot_permutations.R

This needs explaining.


