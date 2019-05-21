context("fit_non_linear")

test_that("fit works", {
  
  test_annotation = example_sample_annotation[example_sample_annotation$MS_batch == 'Batch_1', ]
  selected_files = test_annotation$FullRunName
  
  df_selected = example_proteome[example_proteome$peptide_group_label == example_proteome$peptide_group_label[1],]
  df_selected = df_selected[df_selected$FullRunName %in% selected_files, ]
  df_selected = merge(df_selected, test_annotation, by =  'FullRunName')

  fit_values = fit_nonlinear(df_selected)
  
  expect_length(fit_values, nrow(df_selected))
  #TODO: write additional tests
})
