# Workflow based on example proteome, sample anotation and peptide annotation files in proBatch package 

#### install dependencies ####
#library(ProjectTemplate)
#load.project()

source("http://bioconductor.org/biocLite.R") # from README.md
bioc_deps <- c("GO.db", "preprocessCore", "impute", "sva", "pvca")
cran_deps <- c("tidyverse", "reshape2", "lazyeval", "readr", "WGCNA", "rlang", "corrplot", "ggfortify", "pheatmap")
biocLite(bioc_deps); 
install.packages(cran_deps); 
install.packages("roxygen2"); 
install.packages("lubridate"); 

lapply(bioc_deps, require, character.only = TRUE)
lapply(cran_deps, require, character.only = TRUE)
require(roxygen2)
require(lubridate) 


#### Load data 
load("~/R/proBatch/data/example_peptide_annotation.RData")
load("~/R/proBatch/data/example_sample_annotation1.rda")
load("~/R/proBatch/data/example_proteome.RData")

load("//pasteur/SysBC-Home/leech/Desktop/proBatch/data/example_peptide_annotation.RData")
load("//pasteur/SysBC-Home/leech/Desktop/proBatch/data/example_sample_annotation.rda")
load("//pasteur/SysBC-Home/leech/Desktop/proBatch/data/example_proteome.RData")

sample_annotation = example_sample_annotation1
peptide_annotation = example_peptide_annotation
proteome_dataset = example_proteome

# generate peptide annotation file from the proteome data 
peptide_annotation = create_peptide_annotation(proteome_dataset, 
                                               peptide_col = "peptide_group_label", 
                                               protein_col = c("Uniprot_ID", "Gene")) 

# remove peptides with missing batch - works 
proteome_dataset = remove_peptides_with_missing_batch(proteome_dataset, sample_annotation)

# clean requants - DEBUG!! this changes the fullrunname for some reason!!! 
proteome_dataset = clean_requants(proteome_dataset, sample_annotation, peptide_annotation,
                                  batch_col = 'MS_batch.final',
                                  feature_id_col = 'peptide_group_label',
                                  missing_frac_batch = .3, missing_frac_total = .3)

# convert proteome into matrix - works 
data_matrix = convert_to_matrix(proteome_dataset, feature_id_col = 'peptide_group_label',
                                measure_col = 'Intensity',
                                sample_id_col = 'FullRunName')

# convert the matrix into long data(if necessary) - works 
data_long = matrix_to_long(data_matrix, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')

# make order in sample_annotation file - works (fixed)
sample_annotation = date_to_sample_order (sample_annotation,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          order_col = 'order',
                                          instrument_col = NULL)


############# Data transformation, normalization and batch correction ##############
sample_annotation = fix_sample_annotation
data_matrix = fix_data_matrix

# log2 transformation of data_matrix 
data_matrix_log2 = log2(data_matrix + 1) 

# quantile normalization of the log2 transformed data 
data_matrix_qnorm = quantile_normalize(data_matrix_log2)

# Normalize with LOESS fitting 
log2qnormfit = normalize_custom_fit(data_matrix_qnorm, sample_annotation,
                                    batch_col = 'MS_batch.final',
                                    feature_id_col = 'peptide_group_label',
                                    sample_id_col = 'FullRunName',
                                    measure_col = 'Intensity',
                                    sample_order_col = 'order',
                                    fit_func = fit_nonlinear,
                                    fitFunc = 'loess_regression')
data_matrix_fit = log2qnormfit$data_matrix
data_long_fit = matrix_to_long(data_matrix_fit, feature_id_col = 'peptide_group_label',
                               measure_col = 'Intensity',
                               sample_id_col = 'FullRunName')


# Batch correction with LOESS + ComBat 
data_matrix_combat = correct_with_ComBat (data_matrix_fit, sample_annotation,
                                          batch_col = 'MS_batch.final', par.prior = TRUE)

# Median centering
data_long_medianCentering = normalize_medians_batch(data_long_fit, sample_annotation,
                                                    sample_id_col = 'FullRunName',
                                                    batch_col = 'MS_batch.final',
                                                    feature_id_col = 'peptide_group_label',
                                                    measure_col = 'Intensity')
data_matrix_medianCentering = convert_to_matrix(data_long_medianCentering, feature_id_col = 'peptide_group_label',
                                                measure_col = 'Intensity',
                                                sample_id_col = 'FullRunName')





########### Plot diagnostics for each of the steps ###################
data_1 = data_matrix_log2
data_2 = data_matrix_qnorm
data_3 = data_matrix_fit
data_4 = data_matrix_combat

## the matrix has been converted back to df_long for some diagnostics that require long data frames 
df_long_1 = matrix_to_long(data_1, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName' )
df_long_2 = matrix_to_long(data_2, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')
df_long_3 = matrix_to_long(data_3, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')
df_long_4 = matrix_to_long(data_4, feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', sample_id_col = 'FullRunName')

# plot sample mean - works
plot_sample_mean(data_1, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL, ylimits = c(12,17))

plot_sample_mean(data_2, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

plot_sample_mean(data_3, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)

plot_sample_mean(data_4, sample_annotation = sample_annotation, sample_id_col = 'FullRunName',
                 order_col = 'order',batch_col = "MS_batch.final", facet_col = NULL)


# plot boxplots - works
plot_boxplot(df_long_1,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(df_long_2,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(df_long_3,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)

plot_boxplot(df_long_4,  sample_annotation = sample_annotation,
             sample_id_col = 'FullRunName',
             measure_col = 'Intensity',
             order_col = "order",
             batch_col = 'MS_batch.final',
             facet_col = NULL,
             color_by_batch = T, color_scheme = 'brewer',
             theme = 'classic',
             plot_title = NULL, order_per_facet = F)


# plot clustering - works
colors_list = sample_annotation_to_colors(sample_annotation,  columns_for_plotting = NULL,
                                          sample_id_col = 'FullRunName',
                                          factor_columns = c('MS_batch.final','EarTag', "Strain", "Diet", "Sex"),
                                          not_factor_columns = c("Age_Days", "order"),
                                          rare_categories_to_other = T,
                                          numerics_to_log = F,
                                          numeric_palette_type = 'brewer',
                                          granularity = 10)

color_alldf = colors_list$color_df
color_df = data.frame(MS_batch = color_df$MS_batch.final,
                      Diet = color_df$Diet,
                      Order = color_df$order)

plot_sample_clustering(data_1, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)
plot_sample_clustering(data_2, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)
plot_sample_clustering(data_3, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)
plot_sample_clustering(data_4, color_df,  distance = "euclidean", agglomeration = 'complete',  
                       label_samples = T, label_font = .2, plot_title = NULL)

# plot heatmap - works with sample_annotation file 
## Check: Adding sample annotation file makes an error: 
##  Error in seq.int(rx[1L], rx[2L], length.out = nb) : 'from' must be a finite number 

plot_heatmap(data_1, sample_annotation = NULL, fill_the_missing = T,
             cluster_rows = T, cluster_cols = F,
             annotation_color_list = colors_list,
             heatmap_color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
             color_for_missing = 'black',
             filename = NA, plot_title = NA)


# plot PCA - works
plot_PCA(data_1, sample_annotation, 
         feature_id_col = 'peptide_group_label',
         color_by = 'Diet',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 

plot_PCA(data_2, sample_annotation, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 

plot_PCA(data_3, sample_annotation, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 

plot_PCA(data_4, sample_annotation, 
         feature_id_col = 'peptide_group_label',
         color_by = 'MS_batch.final',
         PC_to_plot = c(1,2), fill_the_missing = 0,
         colors_for_factor = NULL,
         theme = 'classic',
         plot_title = NULL) 


# plot PVCA - works
pvca1 = plot_pvca(data_1, sample_annotation,
                  sample_id_col = 'FullRunName',
                  feature_id_col = 'peptide_group_label',
                  technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
                  biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
                  fill_the_missing = 0,
                  threshold_pca = .6, threshold_var = .01,
                  colors_for_bars = NULL,
                  theme = 'classic', plot_title = NULL)

pvca2 = plot_pvca(data_2, sample_annotation,
                  sample_id_col = 'FullRunName',
                  feature_id_col = 'peptide_group_label',
                  technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
                  biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
                  fill_the_missing = 0,
                  threshold_pca = .6, threshold_var = .01,
                  colors_for_bars = NULL,
                  theme = 'classic', plot_title = NULL)

pvca3 = plot_pvca(data_3, sample_annotation,
                  sample_id_col = 'FullRunName',
                  feature_id_col = 'peptide_group_label',
                  technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
                  biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
                  fill_the_missing = 0,
                  threshold_pca = .6, threshold_var = .01,
                  colors_for_bars = NULL,
                  theme = 'classic', plot_title = NULL)

pvca4 = plot_pvca(data_4, sample_annotation,
                  sample_id_col = 'FullRunName',
                  feature_id_col = 'peptide_group_label',
                  technical_covariates = c('MS_batch.final', 'ProteinPrepDate'),
                  biological_covariates = c('EarTag','Strain', "Diet", "Sex", "Age_Days"),
                  fill_the_missing = 0,
                  threshold_pca = .6, threshold_var = .01,
                  colors_for_bars = NULL,
                  theme = 'classic', plot_title = NULL)

# plot peptide levels of one protein - works 
plot_peptides_of_one_protein (proteinName = "Haao",  protein_col = "Gene", df_long_1, 
                              sample_annotation,example_peptide_annotation,
                              order_col = 'order',
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              measure_col = 'Intensity',
                              feature_id_col = 'peptide_group_label')


# plot peptide trend - Works
rows = which(peptide_annotation$Gene == "Haao")
peptide_Haao = example_peptide_annotation$peptide_group_label[rows]

plot_peptide_trend (pep_name = peptide_Haao[1], 
                    df_long = df_long_1, 
                    sample_annotation = sample_annotation,
                    order_col = "order",
                    sample_id_col = 'FullRunName',
                    batch_col = 'MS_batch.final',
                    measure_col = 'Intensity',
                    feature_id_col = 'peptide_group_label',
                    geom = c('point', 'line'),
                    color_by_batch = F, facet_by_batch = F,
                    requant = NULL,
                    plot_title = NULL,
                    vline_color ='red',
                    theme = 'classic')


# plot spike-ins - Works 
plot_spike_ins_trend(df_long = df_long_1,
                     sample_annotation = sample_annotation,
                     peptide_annotation = peptide_annotation,
                     protein_col = 'Gene',
                     order_col = 'order',
                     spike_ins = "BOVINE_A1ag" ,
                     sample_id_col = 'FullRunName',
                     batch_col = 'MS_batch.final',
                     measure_col = 'Intensity',
                     feature_id_col = 'peptide_group_label',
                     requant = NULL,
                     plot_title = 'Spike-in BOVINE protein peptides')

plot_spike_ins_trend(df_long = df_long_1,
                     sample_annotation = sample_annotation,
                     peptide_annotation = example_peptide_annotation,
                     protein_col = 'Gene',
                     order_col = 'order',
                     spike_ins = "BOVINE_FetuinB" ,
                     sample_id_col = 'FullRunName',
                     batch_col = 'MS_batch.final',
                     measure_col = 'Intensity',
                     feature_id_col = 'peptide_group_label',
                     requant = NULL,
                     plot_title = 'Spike-in BOVINE protein peptides')

#### workflow to generate corrletion-based plots ####
# Consider of having 3 options: corrplot, corrplot.mixed and pheatmap? 
# vignette: mention possible options for each argument e.g. flavor - pheatmap vs. corrplot 
# vignette: illustrate different options for pheatmap and corrplot 

plot_protein_corrplot(data_matrix = data_1,
                      protein_name = "Haao",
                      peptide_annotation = peptide_annotation,
                      protein_col = 'Gene',
                      peptide_col_name = 'peptide_group_label',
                      flavor = 'pheatmap',
                      filename = NULL,
                      width = NA, height = NA, unit = c('cm','in','mm'),
                      plot_title = 'peptide correlation matrix')

plot_protein_corrplot(data_matrix = data_3,
                      protein_name = "Haao",
                      peptide_annotation = peptide_annotation,
                      protein_col = 'Gene',
                      peptide_col_name = 'peptide_group_label',
                      flavor = 'corrplot', 
                      lower.col = "black", number.cex = .5,tl.cex = .7,
                      filename = NULL,
                      width = NA, height = NA, unit = c('cm','in','mm'),
                      plot_title = 'peptide correlation matrix')

#	Plot_samples_corrplot - works 
plot_samples_corrplot(data_matrix = data_1, 
                      samples_to_plot = colnames(data_1)[1:50] ,
                      flavor = 'pheatmap', filename = NULL,
                      width = NA, height = NA, unit = c('cm','in','mm'),
                      plot_title = 'Correlation matrix of samples')

plot_samples_corrplot(data_matrix = data_1, 
                      samples_to_plot = colnames(data_1)[1:10] ,
                      flavor = 'corrplot', tl.cex = .8,
                      lower.col = "black", number.cex = .7, upper = "square", 
                      filename = NULL,
                      width = NA, height = NA, unit = c('cm','in','mm'),
                      plot_title = 'Correlation matrix of samples')

#	Plot_sample_corr_distribution - works
plot_sample_corr_distribution(data_matrix = data_1, sample_annotation,
                              repeated_samples = NULL,
                              sample_id_col = 'FullRunName',
                              batch_col = 'MS_batch.final',
                              biospecimen_id_col = 'EarTag',
                              plot_title = 'Correlation distribution',
                              plot_param = 'batch_replicate')

#	Plot_prot_corr_distribution - works 
plot_prot_corr_distribution(data_matrix = data_1, peptide_annotation,
                            protein_col = 'Gene',
                            feature_id_col = 'peptide_group_label',
                            plot_title = 'Distribution of peptide correlation',
                            theme = 'classic')

#	Plot_within_prot_corr_distribution - works
# Comment: very computationally expensive
matrice_list = join_data_matrices (matrices_list = c(log = data_1, combat = data_4), step = NULL,
                                   sample_annotation, feature_id_col = 'peptide_group_label',
                                   measure_col = 'Intensity',sample_id_col = 'FullRunName')

matrice_list = list(log = data_1, combat = data_4)
plot_within_prot_corr_distribution(matrice_list, peptide_annotation,
                                   protein_col = 'Gene',
                                   feature_id_col = 'peptide_group_label',
                                   plot_title = 'Distribution of peptide correlation',
                                   theme = 'classic')
# WARNING: correlation distribution looks worse after normalization. Need to double check! 

#	Plot_peptide_correlation_distr_one_protein - works 
plot_peptide_correlation_distr_one_protein(matrice_list,
                                           protein_name = "Haao",
                                           peptide_annotation,
                                           protein_col = 'Gene',
                                           feature_id_col = 'peptide_group_label',
                                           plot_title = "Distribution of peptide correlation",
                                           theme = 'classic')




