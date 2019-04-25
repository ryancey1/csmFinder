<div align=center><img width="300" height="320" src="https://github.com/Gavin-Yinld/csmFInder/blob/master/figures/csmFinder.gif"/></div>

# csmFinder

# Introduction

csmFinder is an R package for identifying putative cell-subset specific methylation (pCSM) loci from methylation datasets generated by single cells or bulk tissue. For single cell methylomes, it uses beta mixture model to identify the genomic loci with bipolar methylation pattern across single cells. For bulk methylomes, a nonparametric Bayesian clustering algorithm is used for grouping the sequence reads into hyper- and hypo-methylated subset and determining the genomic loci with significant difference bwtween two subsets as so called pCSM loci. 

The package includes two main function modules, the first one identify pCSM loci from bismark output file. The other perform co-methylation analysis and extract eigen-pCSM loci from each co-methylation module.

# Current Features
* Generate 4-CpG segments from bismark extractor results
* Identify candidate segments covered by totally methylated and unmethylated reads (or single cells)
* Identify pCSM 4-CpG segments with bipolar methylation pattern
* Merge pCSM segments to pCSM loci
* Calculate methylation level in pCSM loci
* Co-methylation analysis to cluster pCSM loci with similar methylation pattern into co-methylated modules
* PCA analysis to extract eigen-pCSM loci representing methylation trend of ecah co-methylation module

# Installation
csmFinder needs the following tools to be installed and available in the `PATH` environment:
1.  [python2](https://www.python.org/downloads/) (>=2.7.10), to process the bismark extractor results
2.  [bedtools2](https://github.com/arq5x/bedtools2) (>=2.27.1) to merge the overlapped pCSM segments into pCSM loci

In R console,
```R
library("devtools")
install_github("Gavin-Yinld/csmFinder")
```
