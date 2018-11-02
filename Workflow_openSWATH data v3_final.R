# FullDataSet Finalized normalization workflow for Evan's analysis #

##################### Install dependencies if necessary ###########################
bioc_deps <- c("GO.db", "preprocessCore", "impute", "sva", "pvca")
cran_deps <- c("roxygen2","lubridate","tidyverse","reshape2","lazyeval","readr","WGCNA", "rlang", 
               "corrplot","ggfortify","pheatmap","dplyr","data.table","wesanderson")
source("http://bioconductor.org/biocLite.R") 
biocLite(bioc_deps); 
install.packages(cran_deps); 
lapply(bioc_deps, require, character.only = TRUE)
lapply(cran_deps, require, character.only = TRUE)

sapply(.libPaths(), function(elt) "AnnotationDbi" %in% dir(elt))
biocLite("GO.db", INSTALL_opts = "--no-test-load")
BiocManager::install("GO.db")

###################################################################################

# Loading data
Fullruns_dataset = as.data.frame(fread(file = 'S:/E1801171630_feature_alignment_requant.tsv', sep = '\t', header = TRUE))
Fullruns_dataset = Fullruns_dataset %>% 
  mutate(FullRunName = gsub('/scratch/55808263.tmpdir/wevan_(.+)\\.mzXML\\.gz', '\\1', filename))

# Peptide annotation  
peptide_annotation = read.csv("S:/html/peptide_annotations_6600evosep.csv")

# Sample annotation 
sample_annotation = read.csv("S:/html/sample_annotation_6600eksigent_REORDER.csv")
sample_annotation = date_to_sample_order (sample_annotation,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          order_col = 'order',
                                          instrument_col = NULL)

# Generate data_matrix from SWATH_long
SWATH_matrix = long_to_matrix(Fullruns_dataset, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity',
                                 sample_id_col = 'FullRunName')

############################### Normalization & Batch correction ####################################################
# Quantile normalize log2 transformed matrix 
normalized_matrix = proBatch::normalize(SWATH_matrix, normalizeFunc = "quantile", log = 2)

# Batch effect correction 
batch_corrected = correct_batch_trend(data_matrix = normalized_matrix, sample_annotation, fitFunc = 'loess_regression', 
                                discreteFunc = 'MedianCentering', batch_col = 'MS_batch.final',  
                                feature_id_col = 'peptide_group_label', sample_id_col = 'FullRunName',
                                measure_col = 'Intensity',  sample_order_col = 'order', 
                                loess.span = 0.75, abs.threshold = 5, pct.threshold = 0.20)

############################# Preliminary plots to confirm the workflow #############################################
SWATH_long = matrix_to_long(SWATH_matrix, feature_id_col = 'peptide_group_label',
                                      measure_col = 'Intensity', sample_id_col = 'FullRunName' )

batch_corrected_long = matrix_to_long(batch_corrected, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName' )

# diagnostic plots 
plot_peptides_of_one_protein (proteinName = "Haao",  protein_col = "Gene", df_long = SWATH_long, 
                              sample_annotation, peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity',
                              feature_id_col = 'peptide_group_label')

plot_peptides_of_one_protein (proteinName = "Haao",  protein_col = "Gene", df_long = batch_corrected_long, 
                              sample_annotation, peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity',
                              feature_id_col = 'peptide_group_label')


