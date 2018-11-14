# proBatch

## General Overview

The proBatch package contains functions for diagnosing and removing batch effects and other unwanted technical variation from high-thoughput experiments. Although the package has primarily been developed for DIA (SWATH) proteomics data, it should also be applicable to most omic data with minor adaptations.
    
As discussed in the corresponding vignettes, the diagnostic part of the package can be broadly divided in:

1. Proteome-wide and 
2. Feature-level functions

In addition to the diagnostic functions we provide a few convenient wrappers for common batch-effect removal approaches such as ComBat and median centering. Furthermore, the package includes non-linear fitting based approaches to deal with complex, mass spectrometry-specific signal drifts.

## Installing

Install dependencies:

```
source("http://bioconductor.org/biocLite.R")

bioc_deps <- c("GO.db", "impute", "preprocessCore", "pvca","sva" )
cran_deps <- c("corrplot", "data.table", "ggfortify","lazyeval", "pheatmap", "reshape2", "rlang", 
                "tidyverse","wesanderson","WGCNA") 

biocLite(bioc_deps) 
install.packages(cran_deps)
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