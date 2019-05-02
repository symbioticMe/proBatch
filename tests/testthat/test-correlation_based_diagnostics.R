context("correlation-based_diagnostics")


test_that("corr_matrix_plots", {
  data(example_proteome_matrix, package="proBatch")
  peptides <- c("10231_QDVDVWLWQQEGSSK_2", "10768_RLESELDGLR_2")
  
  matrix = example_proteome_matrix[peptides,]
  corr_matrix = cor(t(matrix), use = 'complete.obs')
  cor_matrix <- plot_corr_matrix(corr_matrix,  flavor = "corrplot")
  
  expect_equivalent(cor_matrix[1], 1)
  expect_equivalent(cor_matrix[2], 0.9111542373)
  expect_equivalent(cor_matrix[4], 1)
})


test_that("protein_corrplot_plots", {
  data(example_proteome_matrix, package="proBatch")
  data(example_peptide_annotation, package = "proBatch")

  corrplot <- plot_protein_corrplot(example_proteome_matrix, protein_name = 'Haao',
                        peptide_annotation = example_peptide_annotation, 
                        protein_col = 'Gene', flavor = "pheatmap")
  
  expect_equivalent(corrplot$tree_row$method, "complete")
  expect_equivalent(corrplot$tree_row$dist.method, "euclidean")
  
  expect_equivalent(corrplot$tree_row$labels[1], "10231_QDVDVWLWQQEGSSK_2")
  expect_equivalent(corrplot$tree_row$labels[2], "10768_RLESELDGLR_2")
  expect_equivalent(corrplot$tree_row$labels[3],  "1131_AQGSVALSVTQDPAR_2" )

})


test_that("sample_corr_heatmaps", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package = "proBatch")
  
  specified_samples = example_sample_annotation$FullRunName[
    which(example_sample_annotation$order %in% 110:115)] 
  
  expect_warning(sample_heatmap <- plot_sample_corr_heatmap(example_proteome_matrix, 
                           samples_to_plot = specified_samples, 
                           flavor = 'pheatmap',  cluster_rows= TRUE, cluster_cols=TRUE,
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


test_that("sample_distribution_plots", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package = "proBatch")
  
  matrix <- example_proteome_matrix[1:20, ]
  sample_dist <- plot_sample_corr_distribution(matrix,
                                example_sample_annotation, batch_col = 'MS_batch', 
                                biospecimen_id_col = "EarTag", 
                                plot_param = 'batch_replicate')
  
  expect_equivalent(sample_dist$labels$x, "batch_replicate")
  expect_equivalent(sample_dist$labels$y, "correlation")
  
  expect_equivalent(sample_dist$plot_env$batch_col, "MS_batch")
  expect_equivalent(sample_dist$plot_env$biospecimen_id_col, "EarTag")
  expect_equivalent(sample_dist$plot_env$plot_param,  "batch_replicate")
  expect_equivalent(sample_dist$plot_env$repeated_samples, NULL)

})


test_that("peptide_distribution_plots", {
  data(example_proteome_matrix, package="proBatch")
  data(example_peptide_annotation, package = "proBatch")

  matrix <- example_proteome_matrix[1:20, ]
  peptide_dist <- plot_peptide_corr_distribution(matrix, 
                                      example_peptide_annotation, protein_col = 'Gene')
  
  expect_equivalent(peptide_dist$labels$x, "same_protein")
  expect_equivalent(peptide_dist$labels$y, "correlation")
  
  expect_equivalent(peptide_dist$plot_env$protein_col, "Gene")
  expect_equivalent(peptide_dist$plot_env$feature_id_col, "peptide_group_label")
  
})
