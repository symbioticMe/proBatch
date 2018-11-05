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
  mutate(FullRunName = gsub(".*wevan_","",filename) ) %>%
  mutate(FullRunName = gsub('.mzXML\\.gz', '\\1', FullRunName))

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
SWATH_log2 = log_transform(SWATH_matrix)

############################### Normalization & Batch correction ####################################################
# Quantile normalize log2 transformed matrix 
normalized_matrix = proBatch::normalize(SWATH_matrix, normalizeFunc = "quantile", log = 2)

# Batch effect correction 
batch_corrected = proBatch::correct_batch_trend(data_matrix = normalized_matrix, sample_annotation, fitFunc = 'loess_regression', 
                                discreteFunc = 'MedianCentering', batch_col = 'MS_batch.final',  
                                feature_id_col = 'peptide_group_label', sample_id_col = 'FullRunName',
                                measure_col = 'Intensity',  sample_order_col = 'order', 
                                loess.span = 0.75, abs.threshold = 5, pct.threshold = 0.20)

write.csv(normalized_matrix, "log2_qnorm_FullRuns_matrix.csv")
write.csv(batch_corrected, "qnorm_loess_medianCentered_FullRuns_matrix.csv")
write.csv(batch_corrected_long, "qnorm_loess_medianCentered_FullRuns_long.csv")


############################# Preliminary plots to confirm the workflow #############################################
SWATH_long = matrix_to_long(SWATH_matrix, feature_id_col = 'peptide_group_label',
                                      measure_col = 'Intensity', sample_id_col = 'FullRunName' )

SWATHlog2_long = matrix_to_long(SWATH_log2, feature_id_col = 'peptide_group_label',
                            measure_col = 'Intensity', sample_id_col = 'FullRunName' )

normalized_long = matrix_to_long(normalized_matrix, feature_id_col = 'peptide_group_label',
                            measure_col = 'Intensity', sample_id_col = 'FullRunName' )

batch_corrected_long = matrix_to_long(batch_corrected, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName' )


# diagnostic plots 
plot_sample_mean(SWATH_log2, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL, ylimits = c(12, 16))

plot_sample_mean(batch_corrected, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL, ylimits = c(12,16))


plot_peptides_of_one_protein (proteinName = "BOVINE_A1ag",  protein_col = "Gene", df_long = SWATHlog2_long, 
                              sample_annotation, peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity',
                              feature_id_col = 'peptide_group_label')

plot_peptides_of_one_protein (proteinName = "BOVINE_A1ag",  protein_col = "Gene", df_long = batch_corrected_long, 
                              sample_annotation, peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity',
                              feature_id_col = 'peptide_group_label')


plot_PCA(SWATH_log2, sample_annotation, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 

plot_PCA(batch_corrected, sample_annotation, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 
