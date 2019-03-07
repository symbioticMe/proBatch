context("auxiliary")


test_that("long_conversion_to_matrix", {
  data(example_proteome, package="proBatch")

  rows <- c(which(example_proteome$peptide_group_label == "10062_NVGVSFYADKPEVTQEQK_3"),
            which(example_proteome$peptide_group_label == "101233_QGFNVVVESGAGEASK_2"))
  
  proteome <- example_proteome[rows,]
  proteome_matrix <- long_to_matrix(proteome)
  
  expect_equal(rownames(proteome_matrix), 
               c("10062_NVGVSFYADKPEVTQEQK_3" ,"101233_QGFNVVVESGAGEASK_2" ))
  expect_equal(ncol(proteome_matrix), nrow(proteome)/2)
  
})


test_that("matrix_conversion_to_long", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  rows <- c(which(rownames(example_proteome_matrix) == "10062_NVGVSFYADKPEVTQEQK_3"),
            which(rownames(example_proteome_matrix) == "101233_QGFNVVVESGAGEASK_2"))
  
  peptide_matrix <- as.matrix(example_proteome_matrix[rows,])
  proteome_long <- matrix_to_long(peptide_matrix, example_sample_annotation)
  
  expect_equal(unique(proteome_long$peptide_group_label), 
               c("10062_NVGVSFYADKPEVTQEQK_3" ,"101233_QGFNVVVESGAGEASK_2" ))
  expect_equal(nrow(proteome_long), ncol(peptide_matrix)*2)

  setequal(colnames(proteome_long), 
               c(colnames(example_sample_annotation), "peptide_group_label", "Intensity"))
  
})


test_that("generated_peptide_annotation", {
  data(example_proteome, package="proBatch")

  peptide_annt <- create_peptide_annotation(example_proteome, 
                                            feature_id_col = "peptide_group_label", 
                                            annotation_col = c("ProteinName"))
  
  expect_equal(colnames(peptide_annt), c( "peptide_group_label", "ProteinName" ))

})


