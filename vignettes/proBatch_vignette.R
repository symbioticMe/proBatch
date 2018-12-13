## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE, fig.pos = 'h')

## ----setup, include = FALSE----------------------------------------------
chooseCRANmirror(graphics=FALSE, ind=1)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----batch_workflow, include = TRUE, fig.align = "center", echo=FALSE, fig.cap="proBatch in batch correction workflow", out.width = '50%'----
knitr::include_graphics("Batch_effects_workflow_staircase.png")

## ----dependencies, eval = FALSE------------------------------------------
#  bioc_deps <- c("GO.db", "impute", "preprocessCore", "pvca","sva" )
#  cran_deps <- c("corrplot", "data.table", "ggfortify","lazyeval", "pheatmap", "reshape2",
#                 "rlang", "tidyverse","wesanderson","WGCNA")
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(bioc_deps)
#  install.packages(cran_deps)

## ----install_proBatch, fig.show='hold', eval = FALSE---------------------
#  # Once the proBatch package is in Bioconductor, can easily install by:
#  install.packages("proBatch")
#  
#  # Alternatively, install the development version from GitHub:
#  install.packages("devtools")
#  devtools::install_github("symbioticMe/proBatch", build_vignettes = TRUE)

## ----load_tidyverse------------------------------------------------------
require(tidyverse)

## ----col_names-----------------------------------------------------------
feature_id_col = 'peptide_group_label'
measure_col = 'Intensity'
sample_id_col = 'FullRunName'
essential_columns = c(feature_id_col, measure_col, sample_id_col)

## ----tech_bio_cols-------------------------------------------------------
technical_covariates = c('MS_batch', 'digestion_batch', 'RunDate', 'RunTime')
biological_covariates = c('Strain', 'Diet', 'Sex', 'Age_Days')
biospecimen_id_col = "EarTag"

## ----load_data, fig.show='hold'------------------------------------------
library(proBatch)
data("example_proteome", "example_sample_annotation", "example_peptide_annotation", 
     package = "proBatch")

## ----date_to_order, fig.show='hold'--------------------------------------
generated_sample_annotation <- date_to_sample_order(example_sample_annotation,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'generated_DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          new_order_col = 'generated_order',
                                          instrument_col = NULL)
library(knitr)
kable(generated_sample_annotation[1:5,] %>%
  select(c("RunDate", "RunTime", "order", "generated_DateTime", "generated_order")))

## ----pep_annotation, fig.show='hold'-------------------------------------
generated_peptide_annotation <- create_peptide_annotation(example_proteome, 
                                        feature_id_col = 'peptide_group_label',
                                        annotation_col = c('ProteinName', 'Gene'))

## ------------------------------------------------------------------------
example_proteome = example_proteome %>% select(one_of(essential_columns))
gc()

