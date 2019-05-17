context("correct_batch_effects")


test_that("center_feature_batch_medians", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  rows <- which(example_proteome$peptide_group_label == "10062_NVGVSFYADKPEVTQEQK_3")
  
  proteome <- example_proteome[rows,]
  median_proteome <- center_feature_batch_medians(proteome, example_sample_annotation)
  
  n_batch <- length(unique(median_proteome$MS_batch))
  expect_equal(length(unique(median_proteome$diff)), n_batch)
  expect_equal(length(unique(median_proteome$median_batch)), n_batch)
  
})


test_that("adjust_batch_trend", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  short_df <- example_proteome[example_proteome[['peptide_group_label']]
                               %in% c( "10062_NVGVSFYADKPEVTQEQK_3","101233_QGFNVVVESGAGEASK_2"), ]
  
  adjusted <- adjust_batch_trend(short_df, example_sample_annotation, span = 0.7, 
                     abs_threshold = 5, pct_threshold = 0.20)

  n_batch <- length(unique(example_sample_annotation$MS_batch))
  
  expect_equal(adjusted$corrected_df[['peptide_group_label']][1],"10062_NVGVSFYADKPEVTQEQK_3")
  expect_equal(length(unique(adjusted$fit_df$MS_batch)), n_batch)
  expect_equal(adjusted$fit_df$fit[1],2342360.4917360563)

})


test_that("correct_with_ComBat_df", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  short_df <- example_proteome[example_proteome[['peptide_group_label']]
                               %in% c( "10062_NVGVSFYADKPEVTQEQK_3","101233_QGFNVVVESGAGEASK_2"), ]
  combat_df <- correct_with_ComBat_df(short_df, example_sample_annotation)
  
  expect_equal(combat_df[['peptide_group_label']][1],"10062_NVGVSFYADKPEVTQEQK_3")
  expect_equal(combat_df[['Intensity']][1],768661.4)
  
  batch_1 <- example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_1"] 
  batch_2 <- example_sample_annotation$FullRunName[example_sample_annotation$MS_batch == "Batch_2"] 
  
  matrix_batch_1 <- short_df[short_df$FullRunName %in% batch_1, ]
  matrix_batch_2 <- short_df[short_df$FullRunName %in% batch_2, ]
  
  combat_batch_1 <- combat_df[combat_df$FullRunName %in% batch_1, ]
  combat_batch_2 <- combat_df[combat_df$FullRunName %in% batch_2, ]
  
  t_test_matrix <- t.test(matrix_batch_1$Intensity, matrix_batch_2$Intensity)
  t_test_combat <- t.test(combat_batch_1$Intensity, combat_batch_2$Intensity)
  
  expect_lt(t_test_matrix$p.value, 0.05)
  expect_gt(t_test_combat$p.value, 0.05)
    
})
