# FullDataSet Finalized normalization workflow for Evan's analysis #

# Install dependencies if necessary 
bioc_deps <- c("GO.db", "preprocessCore", "impute", "sva", "pvca")
cran_deps <- c("roxygen2","lubridate","tidyverse","reshape2","lazyeval","readr","WGCNA", "rlang", 
               "corrplot","ggfortify","pheatmap","dplyr","data.table","wesanderson")
source("http://bioconductor.org/biocLite.R") 
biocLite(bioc_deps); 
install.packages(cran_deps); 
lapply(bioc_deps, require, character.only = TRUE)
lapply(cran_deps, require, character.only = TRUE)

# Loading data
Fullruns_dataset = as.data.frame(fread(file = 'S:/E1801171630_feature_alignment_requant.tsv', sep = '\t', header = TRUE))
Fullruns_dataset = Fullruns_dataset %>% 
  mutate(FullRunName = gsub('/scratch/55808263.tmpdir/wevan_(.+)\\.mzXML\\.gz', '\\1', filename))

# Peptide annotation  
peptide_annotation = read.csv("S:/html/peptide_annotations_6600evosep.csv")

# Sample annotation 
sample_annotation = read.csv("S:/html/sample_annotation_6600evosep.csv")
sample_annotation = date_to_sample_order (sample_annotation,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          order_col = 'order',
                                          instrument_col = NULL)

# Generate data_matrix from SWATH_long
SWATH_matrix = convert_to_matrix(Fullruns_dataset, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity',
                                 sample_id_col = 'FullRunName')

# For Evan's data, no requant nor peptide filtering are executed. 
# All normalization and batch correction have their intrinsic filtering system for those peptides 
# not suitable for fitting/ normalization/ batch correction. Full dataset 

############################### Normalization & Batch correction ####################################################
# log2 transformation of data_matrix 
SWATH_matrix_log2 = log_transform(SWATH_matrix)

# Quantile normalization of the log2 transformed data 
SWATH_matrix_qnorm = quantile_normalize(SWATH_matrix_log2)

# Batch effect correction 
batch_corrected = correct_batch(data_matrix = SWATH_matrix_qnorm, sample_annotation, fitFunc = 'loess_regression', 
                                discreteFunc = 'ComBat', batch_col = 'MS_batch.final',  
                                feature_id_col = 'peptide_group_label', sample_id_col = 'FullRunName',
                                measure_col = 'Intensity',  sample_order_col = 'order', 
                                loess.span = 0.75, abs.threshold = 5, pct.threshold = 0.30)

############################# Preliminary plots to confirm the workflow #############################################

