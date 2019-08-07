context("correlation-based_diagnostics")


test_that("corr_matrix_plots", {
  data(example_proteome_matrix, package="proBatch")
  peptides <- c("10231_QDVDVWLWQQEGSSK_2", "10768_RLESELDGLR_2")
  
  matrix_test = example_proteome_matrix[peptides,]
  corr_matrix = cor(t(matrix_test), use = 'complete.obs')
  corr_matrix_pheatmap <- plot_corr_matrix(corr_matrix)
  
  expect_is(corr_matrix_pheatmap, 'pheatmap')
  expect_equivalent(corr_matrix_pheatmap$gtable$layout$name[4], "legend")
})


test_that("protein_corrplot_plots", {
  data(example_proteome_matrix, package="proBatch")
  data(example_peptide_annotation, package = "proBatch")

  corrplot <- plot_protein_corrplot(example_proteome_matrix, protein_name = 'Haao',
                        peptide_annotation = example_peptide_annotation, 
                        protein_col = 'Gene', 
                        cluster_rows = T, cluster_cols = T)
  
  expect_equivalent(corrplot$tree_row$method, "complete")
  expect_equivalent(corrplot$tree_row$dist.method, "euclidean")
  
  expect_equivalent(corrplot$tree_row$labels[1], "10231_QDVDVWLWQQEGSSK_2")
  expect_equivalent(corrplot$tree_row$labels[2], "10768_RLESELDGLR_2")
  expect_equivalent(corrplot$tree_row$labels[3],  "1131_AQGSVALSVTQDPAR_2" )

})


test_that("sample_corr_heatmap", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package = "proBatch")
  
  specified_samples = example_sample_annotation$FullRunName[
    which(example_sample_annotation$order %in% 110:115)] 
  
  expect_warning(sample_heatmap <- plot_sample_corr_heatmap(example_proteome_matrix, 
                           samples_to_plot = specified_samples, 
                           cluster_rows= TRUE, cluster_cols=TRUE,
                           annotation_names_col = TRUE, annotation_legend = FALSE, 
                           show_colnames = FALSE))
  
  expect_equivalent(sample_heatmap$tree_row$method, "complete")
  expect_equivalent(sample_heatmap$tree_row$dist.method, "euclidean")
  
  expect_equivalent(sample_heatmap$tree_col$method, "complete")
  expect_equivalent(sample_heatmap$tree_col$dist.method, "euclidean")
            
  expect_equivalent(sample_heatmap$tree_row$labels[1], "I171013_BXD61_CD_ET2145_Run113")
  expect_equivalent(sample_heatmap$tree_row$labels[2], "I171013_BXD70_HF_ET1728_Run114")
  expect_equivalent(sample_heatmap$tree_row$labels[3], "I171016_BXD89_CD_ET2078_Run115")
  
  expect_equivalent(sample_heatmap$gtable$layout[[1]], c(1, 2, 4, 4, 4, 3))
  expect_equivalent(sample_heatmap$gtable$layout$name[[4]], "matrix")
  
})


test_that("sample_distribution_plot", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package = "proBatch")
  
  matrix_test <- example_proteome_matrix[1:20, ]
  sample_dist <- plot_sample_corr_distribution(matrix_test,
                                example_sample_annotation, batch_col = 'MS_batch', 
                                biospecimen_id_col = "EarTag", 
                                plot_param = 'batch_replicate')
  
  expect_equivalent(sample_dist$labels$x, "batch_replicate")
  expect_equivalent(sample_dist$labels$y, "correlation")
  
  expect_is(sample_dist$plot_env$corr_distribution, "data.frame")
  expect_equivalent(sample_dist$plot_env$plot_param,  "batch_replicate")
  expect_is(sample_dist$plot_env$gg, 'ggplot')

})

test_that("calculate_sample_corr_distribution", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package = "proBatch")
  
  matrix_test <- example_proteome_matrix[1:20, ]
  corr_distribution = calculate_sample_corr_distr(data_matrix = matrix_test, 
                                                  repeated_samples = NULL,
                                                  sample_annotation = example_sample_annotation,
                                                  biospecimen_id_col = "EarTag", 
                                                  sample_id_col = 'FullRunName', 
                                                  batch_col = 'MS_batch')
  
  expect_is(corr_distribution, "data.frame")
  
  sample_cols = paste('FullRunName', seq_len(2), sep = "_")
  
  expect_equivalent("batch_replicate" %in% names(corr_distribution), TRUE)
  expect_equivalent("correlation" %in% names(corr_distribution), TRUE)
  expect_equivalent("replicate" %in% names(corr_distribution), TRUE)
  expect_equivalent(all(sample_cols %in% names(corr_distribution)), TRUE)
  
  
})


test_that("peptide_distribution_plots", {
  data(example_proteome_matrix, package="proBatch")
  data(example_peptide_annotation, package = "proBatch")

  matrix_test <- example_proteome_matrix[1:20, ]
  peptide_dist <- plot_peptide_corr_distribution(data_matrix = matrix_test, 
                                                 peptide_annotation = example_peptide_annotation, 
                                                 protein_col = 'Gene')
  
  expect_equivalent(peptide_dist$labels$x, NULL)
  expect_equivalent(peptide_dist$labels$y, "correlation")
  
  expect_is(peptide_dist$plot_env$corr_distribution, "data.frame")
  expect_equal(peptide_dist$plot_env$median_same_prot, 0.8010358, tolerance=1e-6)
  
})
