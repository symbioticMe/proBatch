# proBatch

  <!-- badges: start -->
  [![Lifecycle: maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
  [![BioC status](http://www.bioconductor.org/shields/build/release/bioc/proBatch.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/proBatch)
  [![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

  [![Build Status](https://travis-ci.com/symbioticMe/proBatch.svg?branch=master)](https://travis-ci.com/symbioticMe/proBatch.svg)
  
<!-- badges: end -->

## General Overview

The proBatch package facilitates batch effects analysis and correction high-throughput experiments. 
Although the package has primarily been developed for DIA (SWATH) proteomics data, 
it should also be applicable to most omic data with minor adaptations.
    
The package contains functions for diagnostics (proteome/genome-wide and feature-level), 
correction (normalization and batch effects correction) and quality control.

Diagnostics part of the package features unified color scheme for plotting, 
    that allows to produce publication-quality graphs.

Correction functions are convenient wrappers for common normalization and batch 
effects removal approaches such as quantile normalization and median centering. 
Furthermore, the package includes non-linear fitting based approaches to deal 
with complex, mass spectrometry-specific signal drifts.

Quality control step, mostly based on correlation analysis, allows to assess whether 
the correction improved the quality of the data.

All steps of batch effects analysis and correction are illustrated in the vignette,
    using the subset of real-world large-scale dataset.

Please use following manuscript for citation: 
Diagnostics and correction of batch effects in large-scale proteomic studies: a tutorial. Molecular Systems Biology 17, e10240 (2021).


## Installing

Install the dependencies:

```
bioc_deps <- c("GO.db", "impute", "preprocessCore", "pvca","sva" )
cran_deps <- c("corrplot", "data.table", "ggplot2", "ggfortify","lazyeval", "pheatmap", "reshape2", "rlang", 
               "tibble", "dplyr", "tidyr", "wesanderson","WGCNA") 
               
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(bioc_deps) 
install.packages(cran_deps)
```

NOTE: You might need to also install the following linux packages:
`apt-get install libxml2-dev libz-dev`

Optionally also install:

```
install.packages(c("devtools", "roxygen2", "testthat"))
```


Install proBatch from github:

```
library(devtools)
install_github("symbioticMe/proBatch", build_vignettes = TRUE)
```

## Exploring the package

The complete documentation:
```
help(proBatch)
```

Browse the vignette:
```
browseVignettes('proBatch')
```
