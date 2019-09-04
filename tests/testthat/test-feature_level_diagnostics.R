context("feature-level diagnostics")


test_that("single_feature_plot", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  single_feature <- plot_single_feature(feature_name = "46213_NVGVSFYADKPEVTQEQK_2", 
                                        df_long = example_proteome, 
                                        sample_annotation = example_sample_annotation)

  expect_equal(single_feature$plot_env$feature_name,  "46213_NVGVSFYADKPEVTQEQK_2")
  
  expect_equal(single_feature$labels$x, "order")
  expect_equal(single_feature$labels$y, "Intensity")

  expect_equal(single_feature$plot_env$batch_col, "MS_batch")
  expect_equal(single_feature$plot_env$color_by_batch, FALSE)
  expect_equal(single_feature$plot_env$order_col, "order")
  expect_equal(single_feature$plot_env$vline_color, "red")
  
  
})


test_that("peptides_of_one_protein_plot", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  data(example_peptide_annotation, package="proBatch")
  
  peptides_plot <- plot_peptides_of_one_protein (protein_name = "Haao",
                                peptide_annotation = example_peptide_annotation,
                                protein_col = "Gene", df_long = example_proteome, 
                                sample_annotation = example_sample_annotation, 
                                color_by_batch = TRUE,
                                order_col = 'order', sample_id_col = 'FullRunName', 
                                batch_col = 'MS_batch')  
  
  expect_equal(peptides_plot$plot_env$feature_name[1],  "10231_QDVDVWLWQQEGSSK_2")
  expect_equal(peptides_plot$plot_env$feature_name[2],  "10768_RLESELDGLR_2" )

  expect_equal(peptides_plot$labels$x, "order")
  expect_equal(peptides_plot$labels$y, "Intensity")
  
  expect_equal(peptides_plot$plot_env$batch_col, "MS_batch")
  expect_equal(peptides_plot$plot_env$color_by_batch, TRUE)
  expect_equal(peptides_plot$plot_env$order_col, "order")
  expect_equal(peptides_plot$plot_env$vline_color, "red")
  
})


test_that("spike_in_peptides_plot", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  data(example_peptide_annotation, package="proBatch")
  
  spike_in <- plot_spike_in(spike_ins = "BOVINE_A1ag", 
                            peptide_annotation = example_peptide_annotation,
                            protein_col = 'Gene',  df_long = example_proteome,
                            sample_annotation = example_sample_annotation,
                            plot_title = "Spike-in BOVINE protein peptides")
  
  expect_equal(spike_in$plot_env$feature_name[1],  "10062_NVGVSFYADKPEVTQEQK_3")
  expect_equal(spike_in$plot_env$feature_name[2],  "10063_NVGVSFYADKPEVTQEQKK_3")
  
  expect_equal(spike_in$labels$x, "order")
  expect_equal(spike_in$labels$y, "Intensity")
  
  expect_equal(spike_in$plot_env$batch_col, "MS_batch")
  expect_equal(spike_in$plot_env$color_by_batch, FALSE)
  expect_equal(spike_in$plot_env$order_col, "order")
  expect_equal(spike_in$plot_env$vline_color, "red")
  
})





test_that("iRT_peptides_plot", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  data(example_peptide_annotation, package="proBatch")
  
  iRT <- plot_iRT(irt_pattern = 'iRT', 
                  peptide_annotation = example_peptide_annotation,
                  protein_col = 'Gene',
                  df_long = example_proteome,
                  sample_annotation = example_sample_annotation)
  
  expect_equal(iRT$plot_env$feature_name[1],   "1146_ADVTPADFSEWSK_3")
  expect_equal(iRT$plot_env$feature_name[2], "12476_TPVISGGPYEYR_2")
  
  expect_equal(iRT$labels$x, "order")
  expect_equal(iRT$labels$y, "Intensity")
  
  expect_equal(iRT$plot_env$batch_col, "MS_batch")
  expect_equal(iRT$plot_env$color_by_batch, FALSE)
  expect_equal(iRT$plot_env$order_col, "order")
  expect_equal(iRT$plot_env$vline_color, "red")
  
})



test_that("fitting_trend_plots", {
  data(example_proteome, package="proBatch")
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  short_df <- example_proteome[example_proteome$peptide_group_label %in% 
                                 unique(example_proteome$peptide_group_label)[1:3], ]
  loess_fit <- adjust_batch_trend_df(short_df, example_sample_annotation, span = 0.7)
  
  fit_plot <- plot_with_fitting_curve(feature_name = "10062_NVGVSFYADKPEVTQEQK_3", 
                                      fit_df = loess_fit, fit_value_col = 'fit',
                                      df_long = example_proteome, 
                                      sample_annotation = example_sample_annotation)
  
  expect_equal(fit_plot$plot_env$feature_name[1],   "10062_NVGVSFYADKPEVTQEQK_3")

  expect_equal(fit_plot$labels$x, "order")
  expect_equal(fit_plot$labels$y, "Intensity")
  
  expect_equal(fit_plot$plot_env$batch_col, "MS_batch")
  expect_equal(fit_plot$plot_env$color_by_batch, FALSE)
  expect_equal(fit_plot$plot_env$order_col, "order")
  expect_equal(fit_plot$plot_env$vline_color, "grey")
  
})
