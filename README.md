# proBatch

## General Overview

The proBatch package contains functions for diagnosing and removing batch effects and other unwanted sources of variation from high-thoughput experiments. Altough the package hass primarily been developed with DIA proteomics data in mind, it should also be applicable to, similar, type of data.
    
As discussed in the corresponding vignettes, the diagnostic part of the package can be broadly divided in:

1. Proteome-wide and 
2. Protein-specific functions

In addition to the diagnostic functions we provide a few convenience wrappers for common batch-effect removal approaches like ComBat [^1] and mean/median centering. Furthermore, the package includes non-linear fitting based approaches to deal with complex, proteomics-specific technical artifacts.

## Installing

Install dependencies:

```
source("http://bioconductor.org/biocLite.R")

bioc_deps <- c("GO.db", "preprocessCore", "impute", "sva", "pvca")
cran_deps <- c("tidyverse", "reshape2", "lazyeval", "readr", "WGCNA", "rlang", "corrplot", "ggfortify", "pheatmap", "wesanderson")

biocLite(biov_deps) 
install.packages(dep)
```

NOTE: You might need to also install the following linux packages:
`apt-get install libxml2-dev libz-dev`

Optionally also install:

```
install.packages("roxygen2")
```


Install proBatch from github:

```
library(devtools)
install_github("symbioticMe/proBatch.git"")
```


## Getting started

POINT TO THE VIGNETTE
ADD A FEW EXAMPLES

## Citing

To reference the pipeline please cite [^1].
If you use tools of this package in your analysis please reference [^2]. 
    
[^1]: Čuklina et al. 2019, MSB
[^2]: Čuklina et al., 2019, Bioinformatics