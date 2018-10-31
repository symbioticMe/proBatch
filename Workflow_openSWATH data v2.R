# Set library directory  
.libPaths("S:/proBatch/%HOMESHARE%/R3UserLibs")
.libPaths("C:/Users/leech/Documents/proBatch/%HOMESHARE%/R3UserLibs")
.libPaths("\\\\pasteur/sysBC-Home/leech/Documents/proBatch/%HOMESHARE%/R3UserLibs")
setwd("\\\\pasteur/sysBC-Home/leech/Documents/proBatch")

# Install dependencies 
bioc_deps <- c("GO.db", "preprocessCore", "impute", "sva", "pvca")
cran_deps <- c("roxygen2","lubridate","tidyverse","reshape2","lazyeval","readr","WGCNA", "rlang", 
              "corrplot","ggfortify","pheatmap","dplyr","data.table","wesanderson")

source("http://bioconductor.org/biocLite.R") 
biocLite(bioc_deps); 
install.packages(cran_deps); 
lapply(bioc_deps, require, character.only = TRUE)
lapply(cran_deps, require, character.only = TRUE)

# OpenSWATH Data (load data of interest)
Allruns_requant = read.table(file = "S:/html/AllRuns_Evosep_feature_alignment_requant.tsv", sep = '\t', header = TRUE)
Allruns_requant = Allruns_requant %>% mutate(FullRunName = gsub('wevan_(.+)\\.mzXML\\.gz', '\\1', filename))

library(readr) #this is full Evan's dataset
Fullruns_dataset = read_delim("S:/html/E1801171630_feature_alignment_requant.tsv", "\t", escape_double = FALSE)
Fullruns_dataset = Fullruns_dataset %>% mutate(FullRunName = gsub('/scratch/55808263.tmpdir/wevan_(.+)\\.mzXML\\.gz', '\\1', filename))

#Fullrun
Fullruns_dataset = read.table("S:/E1801171630_feature_alignment_requant.tsv", "\t",header = TRUE)
Fullruns_dataset = read.table(file = "S:/html/E1801171630_feature_alignment_requant.tsv", sep = '\t', header = TRUE)

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

############## pre-processing #############################
openSWATH_data = Fullruns_dataset

# remove peptides with missing batch - crashes
openSWATH_data = remove_peptides_with_missing_batch(openSWATH_data, sample_annotation,
                                                            batch_col = 'MS_batch.final',
                                                            feature_id_col = 'peptide_group_label')

############## normalization ###############################
# Generate data_matrix from SWATH_long
SWATH_matrix = convert_to_matrix(Allruns_requant, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity',
                                 sample_id_col = 'FullRunName')

# log2 transformation of data_matrix 
SWATH_matrix_log2 = log2(SWATH_matrix + 1) 

# quantile normalization of the log2 transformed data 
SWATH_matrix_qnorm = quantile_normalize(SWATH_matrix_log2)

# LOESS fitting 
SWATH_list_fit = normalize_custom_fit(SWATH_matrix_qnorm, sample_annotation,
                                      batch_col = 'MS_batch.final',
                                      feature_id_col = 'peptide_group_label',
                                      sample_id_col = 'FullRunName',
                                      measure_col = 'Intensity',
                                      sample_order_col = 'order',
                                      fit_func = fit_nonlinear,
                                      fitFunc = 'loess_regression')
SWATH_matrix_fit = SWATH_list_fit$data_matrix
SWATH_long_fit = matrix_to_long(SWATH_matrix_fit, feature_id_col = 'peptide_group_label',
                                measure_col = 'Intensity', sample_id_col = 'FullRunName')

# Median centering
data_long_medianCentering = normalize_medians_batch(SWATH_long_fit, sample_annotation,
                                                    sample_id_col = 'FullRunName',
                                                    batch_col = 'MS_batch.final',
                                                    feature_id_col = 'peptide_group_label',
                                                    measure_col = 'Intensity')
data_matrix_medianCentering = convert_to_matrix(data_long_medianCentering, feature_id_col = 'peptide_group_label',
                                                measure_col = 'Intensity_normalized',
                                                sample_id_col = 'FullRunName')


# ComBat noramlization (on hold, median centering should be good enough)
SWATH_matrix_ComBat = correct_with_ComBat (SWATH_matrix_fit, sample_annotation,
                                                      batch_col = 'MS_batch.final', par.prior = TRUE)


# prepare full dataset for Wenguang's pipeline as openSWATH output - qnorm + loess + medianCentering 
normalized_long = data.frame(FullRunName = data_long_medianCentering$FullRunName,
                             peptide_group_label = data_long_medianCentering$peptide_group_label,
                             Intensity_normalized = data_long_medianCentering$Intensity_normalized)

noramlized_openSWATH = merge(normalized_long, Allruns_requant)
write.csv(noramlized_openSWATH, "~/R/AllRuns; qnorm_loess_medianCenter.csv", row.names=T)
write.csv(peptide_annotation, "~/R/peptide_annotation_6600.csv", row.names=T)
write.csv(sample_annotation, "~/R/sample_annotation_6600.csv", row.names=T)


################### convert to data_long for plotting #########################
# matrix: SWATH_matrix_log2, SWATH_matrix_qnorm, SWATH_matrix_fit, SWATH_matrix_ComBat

SWATH_long_log2 = matrix_to_long(SWATH_matrix_log2, feature_id_col = 'peptide_group_label',
                            measure_col = 'Intensity', sample_id_col = 'FullRunName')

SWATH_long_qnorm = matrix_to_long(SWATH_matrix_qnorm, feature_id_col = 'peptide_group_label',
                                 measure_col = 'Intensity', sample_id_col = 'FullRunName')

#SWATH_long_fit = matrix_to_long(SWATH_matrix_fit, feature_id_col = 'peptide_group_label',
#                                 measure_col = 'Intensity', sample_id_col = 'FullRunName')

#SWATH_long_ComBat = matrix_to_long(SWATH_matrix_ComBat, feature_id_col = 'peptide_group_label',
#                                 measure_col = 'Intensity', sample_id_col = 'FullRunName')

######################### Data preparation #############################################
df_long_1 = matrix_to_long(data_1, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName' )
df_long_2 = matrix_to_long(data_2, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')
df_long_3 = matrix_to_long(data_3, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')
df_long_4 = matrix_to_long(data_4, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')


################## plotting diagnostics of normalization ###########################

# plot sample mean - works, double check of the change 
plot_sample_mean(SWATH_matrix_log2, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL, ylimits = c(15.5, 17.0))

plot_sample_mean(SWATH_matrix_qnorm, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

plot_sample_mean(SWATH_matrix_fit, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

plot_sample_mean(data_matrix_medianCentering, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

# boxplots 
plot_boxplot(SWATH_long_log2,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(SWATH_long_qnorm,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(SWATH_long_fit,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(data_long_medianCentering,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

# plot peptides of one protein 
plot_peptides_of_one_protein (proteinName = "Haao",  protein_col = "Gene", df_long = SWATH_long_qnorm, 
                              sample_annotation, peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity',
                              feature_id_col = 'peptide_group_label')

plot_peptides_of_one_protein (proteinName = "Haao",  protein_col = "Gene", df_long = data_long_medianCentering, 
                              sample_annotation, peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity_normalized',
                              feature_id_col = 'peptide_group_label')
