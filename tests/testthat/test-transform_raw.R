context("transform_raw")

test_that("log_transformed_matrix", {
  data(example_proteome_matrix, package="proBatch")
  
  matrix_test <- example_proteome_matrix[1:10, ]
  log2_transformed_matrix <- log_transform_matrix(matrix_test, log_base = 2, offset = 0.5)
  log10_transformed_matrix <- log_transform_matrix(matrix_test, log_base = 10, offset = 1)
  
  log2_matrix = log2(matrix_test + 0.5)
  log10_matrix = log10(matrix_test + 1)
  
  expect_equivalent(log2_matrix[1:10], log2_transformed_matrix[1:10])
  expect_equivalent(log10_matrix[1:10], log10_transformed_matrix[1:10])
  
  expect_warning(log_warn <- log_transform_matrix(matrix_test, log_base = NULL))
  
})
