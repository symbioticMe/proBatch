context("correct_batch_effects")


test_that("center_peptide_batch_medians", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  rows <- which(example_proteome$peptide_group_label == "10062_NVGVSFYADKPEVTQEQK_3")
  
  proteome <- example_proteome[rows,]
  median_proteome <- center_peptide_batch_medians(proteome, example_sample_annotation)
  
  n_batch <- length(unique(median_proteome$MS_batch))
  expect_equal(length(unique(median_proteome$diff)), n_batch)
  expect_equal(length(unique(median_proteome$median_batch)), n_batch)
  
})


test_that("adjust_batch_trend", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  rows <- c(which(rownames(example_proteome_matrix) == "10062_NVGVSFYADKPEVTQEQK_3"),
            which(rownames(example_proteome_matrix) == "101233_QGFNVVVESGAGEASK_2"))
  
  peptide_matrix <- as.matrix(example_proteome_matrix[rows,])
  adjusted <- adjust_batch_trend(peptide_matrix, example_sample_annotation, span = 0.7, 
                     abs_threshold = 5, pct_threshold = 0.20)

  n_batch <- length(unique(example_sample_annotation$MS_batch))
  
  expect_equal(rownames(adjusted$data_matrix)[1],"10062_NVGVSFYADKPEVTQEQK_3")
  expect_equal(length(unique(adjusted$fit_df$MS_batch)), n_batch)
  expect_equal(adjusted$fit_df$fit[1],2342360.4917360563)

})


test_that("correct_with_ComBat", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  rows <- c(which(rownames(example_proteome_matrix) == "10062_NVGVSFYADKPEVTQEQK_3"),
            which(rownames(example_proteome_matrix) == "101233_QGFNVVVESGAGEASK_2"))
  
  peptide_matrix <- as.matrix(example_proteome_matrix[rows,])
  combat <- correct_with_ComBat(peptide_matrix, example_sample_annotation)
  
  expect_equal(rownames(combat)[1],"10062_NVGVSFYADKPEVTQEQK_3")
  expect_equal(combat[1,1],768661.4)
  
  batch_1 <- example_sample_annotation$FullRunName[which(example_sample_annotation$MS_batch == "Batch_1")] 
  batch_2 <- example_sample_annotation$FullRunName[which(example_sample_annotation$MS_batch == "Batch_2")] 
  
  matrix_batch_1 <- peptide_matrix[,colnames(peptide_matrix) %in% batch_1]
  matrix_batch_2 <- peptide_matrix[,colnames(peptide_matrix) %in% batch_2]
  
  combat_batch_1 <- combat[,colnames(combat) %in% batch_1]
  combat_batch_2 <- combat[,colnames(combat) %in% batch_2]
  
  t_test_matrix <- t.test(matrix_batch_1,matrix_batch_2)
  t_test_combat <- t.test(combat_batch_1,combat_batch_2)
  
  expect_lt(t_test_matrix$p.value, 0.05)
  expect_gt(t_test_combat$p.value, 0.05)
    
})
