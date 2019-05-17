context("normalize")

test_that("quantile_normalize", {
  data(example_proteome_matrix, package="proBatch")
  
  matrix <- example_proteome_matrix[1:10, ]
  quant_normalized_matrix <- quantile_normalize(matrix)

  quant_matrix <- preprocessCore::normalize.quantiles(matrix)
  
  expect_equivalent(quant_matrix[1:10], quant_normalized_matrix[1:10])
})



test_that("normalize_sample_medians", {
  data(example_proteome, package="proBatch")
  
  sampled_rows <- c( 95576, 51906, 25928, 53092, 12233, 44203, 37086 , 5742, 54461,  6157, 76765, 
             36458, 16956, 83006,  6865, 35359, 34456, 53969, 16449, 34586, 69234, 61117, 25524, 
             54715, 25445, 13009, 71210, 11979, 75665, 48261, 34331, 44354, 31732,  6913, 81975, 
             19950, 14439, 53985, 44994, 8980, 54782, 46463, 67231, 66593, 23441, 50172, 24628, 72037, 49198)
  sample_df <- example_proteome[sampled_rows, ]
  median_centered_data <- normalize_sample_medians(sample_df)

  expect_equivalent(median_centered_data$median_global[1], 32298.40)
  expect_equivalent(median_centered_data$median_global[2], 32298.40)
  
  expect_equivalent(median_centered_data$Intensity[1], 32298.4)
  expect_equivalent(median_centered_data$Intensity[2], 32298.4)
  
})
