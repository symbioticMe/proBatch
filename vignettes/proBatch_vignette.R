## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(echo=TRUE, warning=FALSE, message=FALSE, fig.pos = 'h')

## ----setup, include = FALSE----------------------------------------------
chooseCRANmirror(graphics=FALSE, ind=1)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- include = TRUE, fig.align = "center", echo=FALSE, fig.cap="proBatch in batch correction workflow", out.width = '50%'----
knitr::include_graphics("Batch_effects_workflow_staircase.png")

## ---- eval = FALSE-------------------------------------------------------
#  bioc_deps <- c("GO.db", "impute", "preprocessCore", "pvca","sva" )
#  cran_deps <- c("corrplot", "data.table", "ggfortify","lazyeval", "pheatmap", "reshape2",
#                 "rlang", "tidyverse","wesanderson","WGCNA")
#  source("https://bioconductor.org/biocLite.R")
#  biocLite(bioc_deps)
#  install.packages(cran_deps)

## ---- fig.show='hold', eval = FALSE--------------------------------------
#  # Once the proBatch package is in Bioconductor, can easily install by:
#  install.packages("proBatch")
#  
#  # Alternatively, install the development version from GitHub:
#  install.packages("devtools")
#  devtools::install_github("symbioticMe/proBatch")

## ------------------------------------------------------------------------
require(tidyverse)

## ------------------------------------------------------------------------
feature_id_col = 'peptide_group_label'
measure_col = 'Intensity'
sample_id_col = 'FullRunName'
essential_columns = c(feature_id_col, measure_col, sample_id_col)

## ------------------------------------------------------------------------
technical_covariates = c('MS_batch', 'digestion_batch', 'RunDate', 'RunTime')
biological_covariates = c('Strain', 'Diet', 'Sex', 'Age_Days')
biospecimen_id_col = "EarTag"

## ---- fig.show='hold'----------------------------------------------------
library(proBatch)
data("example_proteome", "example_sample_annotation", "example_peptide_annotation", 
     package = "proBatch")

## ---- fig.show='hold'----------------------------------------------------
generated_sample_annotation <- date_to_sample_order(example_sample_annotation,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'generated_DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          new_order_col = 'generated_order',
                                          instrument_col = NULL)
library(knitr)
kable(generated_sample_annotation[1:5,] %>%
  select(c("RunDate", "RunTime", "order", "generated_DateTime", "generated_order")))

## ---- fig.show='hold'----------------------------------------------------
generated_peptide_annotation <- create_peptide_annotation(example_proteome, 
                                        feature_id_col = 'peptide_group_label',
                                        annotation_col = c('ProteinName', 'Gene'))

## ------------------------------------------------------------------------
example_proteome = example_proteome %>% select(one_of(essential_columns))
gc()

## ------------------------------------------------------------------------
example_matrix <- long_to_matrix(example_proteome)

## ------------------------------------------------------------------------
log_transformed_matrix <- log_transform(example_matrix)

## ---- fig.show='hold'----------------------------------------------------
color_scheme <- sample_annotation_to_colors (example_sample_annotation, 
       factor_columns = c('MS_batch','EarTag', "Strain", "Diet", "digestion_batch", "Sex"),
       not_factor_columns = 'DateTime',
       numeric_columns = c('Age_Days', 'order'))
color_list = color_scheme$list_of_colors

## ---- fig.show='hold', fig.width=5, fig.height=2-------------------------
batch_col = 'MS_batch'
plot_sample_mean(log_transformed_matrix, example_sample_annotation, order_col = 'order', 
                 batch_col = batch_col, color_by_batch = T, ylimits = c(12, 16),
                 color_scheme = color_list[[batch_col]])

## ---- fig.show='hold', fig.width=10, fig.height=5------------------------
log_transformed_long <- matrix_to_long(log_transformed_matrix)
batch_col = 'MS_batch'
plot_boxplot(log_transformed_long, example_sample_annotation, 
             batch_col = batch_col, color_scheme = color_list[[batch_col]])

## ---- fig.show='hold'----------------------------------------------------
median_normalized_matrix = normalize_data(log_transformed_matrix, 
                                          normalizeFunc = "medianCentering")

## ---- fig.show='hold'----------------------------------------------------
median_normalized_matrix = normalize_data(example_matrix, 
                                          normalizeFunc = "medianCentering", log_base = 2)

## ---- fig.show='hold'----------------------------------------------------
quantile_normalized_matrix = normalize_data(log_transformed_matrix, 
                                            normalizeFunc = "quantile")

## ---- fig.show='hold', fig.width=5, fig.height=2-------------------------
plot_sample_mean(quantile_normalized_matrix, example_sample_annotation, 
                 color_by_batch = T, ylimits = c(12, 16), 
                 color_scheme = color_list[[batch_col]])

## ---- fig.show='hold'----------------------------------------------------
batch_corrected_matrix <- correct_batch_effects(data_matrix = quantile_normalized_matrix, 
                                  example_sample_annotation, discreteFunc = 'ComBat',
                                  abs.threshold = 5, pct.threshold = 0.20)

## ---- fig.show='hold', fig.width=10, fig.height=5------------------------
color_annotation <- color_scheme$color_df
selected_annotations <- c("MS_batch",  "digestion_batch", "Diet")
#select only a subset of samples for plotting
color_annotation <- color_annotation %>% select(one_of(selected_annotations))

#Plot clustering between samples 
plot_hierarchical_clustering(quantile_normalized_matrix, color_annotation,  
                             distance = "euclidean", agglomeration = 'complete',
                             label_samples = F)

## ---- fig.show='hold', fig.width=10, fig.height=15-----------------------
plot_heatmap(quantile_normalized_matrix, example_sample_annotation, 
             sample_annotation_col = selected_annotations, 
             cluster_cols = T, 
             annotation_color_list = color_scheme$list_of_colors,
             show_rownames = F, show_colnames = F)

## ---- fig.show='hold',  fig.width=3.4, fig.height=2.3--------------------
plot_PCA(quantile_normalized_matrix, example_sample_annotation, color_by = 'MS_batch', 
              plot_title = "MS batch", colors_for_factor = color_list[[batch_col]])
plot_PCA(quantile_normalized_matrix, example_sample_annotation, color_by = "digestion_batch", 
         plot_title = "Digestion batch", colors_for_factor = color_list[["digestion_batch"]])
plot_PCA(quantile_normalized_matrix, example_sample_annotation, color_by = "Diet",  
         plot_title = "Diet", colors_for_factor = color_list[["Diet"]])

## ---- fig.show='hold', eval = FALSE--------------------------------------
#  pvca <- plot_PVCA(quantile_normalized_matrix, example_sample_annotation,
#                    technical_covariates = c('MS_batch', 'digestion_batch'),
#                    biological_covariates = c(biological_covariates, biospecimen_id_col))

## ---- include = TRUE, fig.align = "center", echo=FALSE, out.width = '80%'----
knitr::include_graphics("pvca_quantile_normalized.png")

## ---- fig.show='hold'----------------------------------------------------
quantile_normalized_long <- matrix_to_long(quantile_normalized_matrix, example_sample_annotation)
plot_spike_in(quantile_normalized_long, example_sample_annotation, 
              peptide_annotation = generated_peptide_annotation,
              protein_col = 'Gene', spike_ins = "BOVINE_A1ag", 
              plot_title = 'Spike-in BOVINE protein peptides',
              color_by_batch = T, color_scheme = color_list[[batch_col]])


## ---- fig.show='hold', fig.width=5, fig.height=2.4-----------------------
loess_fit <- adjust_batch_trend(quantile_normalized_matrix, example_sample_annotation)
loess_fit_matrix <- loess_fit$data_matrix

## ---- fig.show='hold', fig.width=5, fig.height=2.4-----------------------
loess_fit_30 <- adjust_batch_trend(quantile_normalized_matrix, example_sample_annotation, span = 0.3)

quantile_normalized_long <- matrix_to_long(quantile_normalized_matrix)
plot_with_fitting_curve(pep_name = "10231_QDVDVWLWQQEGSSK_2", 
            df_long = quantile_normalized_long, example_sample_annotation, 
            color_by_batch = T, color_scheme = color_list[[batch_col]],
            fit_df = loess_fit_30$fit_df, plot_title = "Span = 30%")

## ---- fig.show='hold', fig.width=5, fig.height=2.4-----------------------
loess_fit_70 <- adjust_batch_trend(quantile_normalized_matrix, example_sample_annotation, span = 0.7)
plot_with_fitting_curve(pep_name = "10231_QDVDVWLWQQEGSSK_2", 
            df_long = quantile_normalized_long, example_sample_annotation, 
            color_by_batch = T, color_scheme = color_list[[batch_col]],
            fit_df = loess_fit_70$fit_df, plot_title = "Span = 70%")

## ------------------------------------------------------------------------
loess_fit_df <- matrix_to_long(loess_fit_matrix)

## ---- fig.show='hold', fig.width=3, fig.height=2.4-----------------------
peptide_median_df <- center_peptide_batch_medians(loess_fit_df, example_sample_annotation)
plot_single_feature(pep_name = "46213_NVGVSFYADKPEVTQEQK_2", df_long = peptide_median_df, 
            example_sample_annotation, color_by_col = NULL, measure_col = 'Intensity_normalized',
            plot_title = "Feature-level Median Centered")

## ---- fig.show='hold'----------------------------------------------------
comBat_matrix <- correct_with_ComBat(loess_fit_matrix, example_sample_annotation)

## ---- fig.show='hold',  fig.width=3, fig.height=2.4----------------------
combat_df <- matrix_to_long(comBat_matrix)
plot_single_feature (pep_name = "46213_NVGVSFYADKPEVTQEQK_2", loess_fit_df, 
          example_sample_annotation, plot_title = "Loess Fitted", color_by_col = NULL)
plot_single_feature (pep_name = "46213_NVGVSFYADKPEVTQEQK_2", combat_df, 
          example_sample_annotation, plot_title = "ComBat corrected", color_by_col = NULL)

## ---- fig.show='hold'----------------------------------------------------
batch_corrected_matrix <- correct_batch_effects(data_matrix = quantile_normalized_matrix, 
                                  example_sample_annotation, discreteFunc = 'ComBat',
                                  abs.threshold = 5, pct.threshold = 0.20)

## ---- fig.show='hold', fig.height=5, fig.width=8-------------------------
earTags <- c("ET1524", "ET2078", "ET1322", "ET1566", "ET1354", "ET1420", "ET2154",
             "ET1515", "ET1506", "ET2577", "ET1681", "ET1585", "ET1518", "ET1906")

# Prepare color annotation 
factors_to_show = c("MS_batch", "EarTag")
replicate_annotation <- example_sample_annotation %>%
  filter(MS_batch == 'Batch_2' | MS_batch == "Batch_3") %>%
  filter(EarTag %in% earTags) %>%
  remove_rownames %>% 
  column_to_rownames(var="FullRunName") %>%
  select(factors_to_show) # Annotate MS_batch and EarTag on pheatmap 

# sample ID of biological replicates 
replicate_filenames = replicate_annotation %>%
  rownames()

breaksList <- seq(0.7, 1, by = 0.01) # color scale of pheatmap 
heatmap_colors = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList))

## ---- fig.show='hold', fig.height=5, fig.width=8-------------------------
# Plot the heatmap 
plot_sample_corr_heatmap(quantile_normalized_matrix, samples_to_plot = replicate_filenames, 
                         flavor = 'pheatmap', plot_title = 'Quantile Normalized', 
                         annotation_colors = color_list[factors_to_show], 
                         annotation_col = replicate_annotation,
                         color = heatmap_colors, breaks = breaksList, 
                         cluster_rows=F, cluster_cols=F,
                         annotation_names_col = TRUE, annotation_legend = FALSE, 
                         show_colnames = F)

## ---- fig.show='hold', fig.height=5, fig.width=8-------------------------
plot_sample_corr_heatmap(batch_corrected_matrix, samples_to_plot = replicate_filenames, 
                         flavor = 'pheatmap', plot_title = 'Batch Corrected',
                         annotation_colors = color_list[factors_to_show], 
                         annotation_col = replicate_annotation,
                         color = heatmap_colors, breaks = breaksList, 
                         cluster_rows=F, cluster_cols=F,
                         annotation_names_col = TRUE, annotation_legend = FALSE, 
                         show_colnames = F)

## ---- fig.show='hold', fig.width=3.2, fig.height=3.5---------------------
sample_cor_norm <- plot_sample_corr_distribution(quantile_normalized_matrix,
                                                 example_sample_annotation, 
                                                 batch_col = 'MS_batch', 
                                                 biospecimen_id_col = "EarTag", 
                                                 plot_title = 'Quantile normalized',
                                                 plot_param = 'batch_replicate')
sample_cor_norm + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0.7,1)

sample_cor_batchCor <- plot_sample_corr_distribution(batch_corrected_matrix,
                                                     example_sample_annotation, 
                                                     batch_col = 'MS_batch', 
                                                     plot_title = 'Batch corrected',
                                                     plot_param = 'batch_replicate')
sample_cor_batchCor + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + ylim(0.7, 1)

## ---- fig.show='hold',fig.width=3, fig.height=2.5------------------------
peptide_cor_norm <- plot_peptide_corr_distribution(quantile_normalized_matrix, 
          generated_peptide_annotation, protein_col = 'Gene', plot_title = 'Quantile normalized')
peptide_cor_norm + geom_hline(yintercept=0, linetype="dashed", color = "grey")

peptide_cor_batchCor <- plot_peptide_corr_distribution(batch_corrected_matrix, 
          generated_peptide_annotation, protein_col = 'Gene', plot_title = 'Batch corrected')
peptide_cor_batchCor + geom_hline(yintercept=0, linetype="dashed", color = "grey")

