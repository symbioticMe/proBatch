context("proteome_wide_diagnostics")


test_that("pca_heatmap", {
  data(example_proteome_matrix, package="proBatch")
  data(sample_color_scheme, package="proBatch")
  
  matrix_test <- example_proteome_matrix[1:10, ]
  color_annotation <- sample_color_scheme$color_df[,c("MS_batch", "Diet")]
 
  hiearchical <- plot_hierarchical_clustering(matrix_test, color_annotation,  
                               distance = "euclidean", agglomeration = 'complete',
                               label_samples = FALSE)
  
  expect_identical(hiearchical$mar[[1]], 1)
  expect_identical(hiearchical$mar[[2]], 5)
  expect_identical(hiearchical$mar[[3]], 0)
  expect_identical(hiearchical$mar[[4]], 1)
  
})


test_that("heatmap_plot", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  data(sample_color_scheme, package="proBatch")
  
  matrix_test <- example_proteome_matrix[1:20, ]
  
  expect_warning(heatmap <- plot_heatmap(matrix_test, sample_annotation = example_sample_annotation, 
                      sample_annotation_col = c("MS_batch",  "digestion_batch", "Diet"), 
                      cluster_cols = TRUE, annotation_color_list = sample_color_scheme$list_of_colors,
                      show_rownames = TRUE, show_colnames = FALSE))
  
  expect_equal(heatmap$tree_row$method, "complete")
  expect_equal(heatmap$tree_row$dist.method , "euclidean")
  
  expect_equal(heatmap$tree_row$labels[1], "10062_NVGVSFYADKPEVTQEQK_3")
  expect_equal(heatmap$tree_row$labels[2], "10063_NVGVSFYADKPEVTQEQKK_3")
  
  expect_equal(heatmap$gtable$layout$name[1], "col_tree")
  expect_equal(heatmap$gtable$layout$name[3], "matrix")
  expect_equal(heatmap$gtable$layout$name[5], "col_annotation")
  expect_equal(heatmap$gtable$layout$name[8], "legend")
  
  expect_equal(heatmap$gtable$layout$t, c(2, 4, 4, 4, 3, 3, 3, 3))

})


test_that("pvca_plot", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  matrix_test <- example_proteome_matrix[1:50, ]
  pvca <- plot_PVCA(matrix_test, example_sample_annotation, 
                   technical_covariates = c('MS_batch', 'digestion_batch'),
                   biological_covariates = c("Diet", "Sex", "Strain"))
  
  expect_equivalent(pvca$plot$data$label[1], factor("Strain"))
  expect_equivalent(pvca$plot$data$label[2], factor("Sex:Strain"))
  expect_equivalent(pvca$plot$data$label[3], factor("digestion_batch"))
  expect_equivalent(pvca$plot$data$label[4], factor("MS_batch"))
  
  expect_equal(pvca$plot$data$category[1], "biological")
  expect_equal(pvca$plot$data$category[3], "technical")
  
})


test_that("pca_plot", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")

  pca <- plot_PCA(example_proteome_matrix, example_sample_annotation, 
                        color_by = 'MS_batch', plot_title = "PCA colored by MS batch")
  expect_equal(pca$labels$y, "PC2 (14.23%)")
  expect_equal(pca$labels$x, "PC1 (69.41%)")
  expect_equal(pca$labels$colour, "MS_batch")
})
