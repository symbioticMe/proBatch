context("initial_assessment")


test_that("sample_mean_plots", {
  data(example_proteome_matrix, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  matrix <- example_proteome_matrix[1:20, ]
  meanplot <- plot_sample_mean(matrix, example_sample_annotation, 
                   order_col = 'order', batch_col = "MS_batch", color_by_batch = TRUE)

  expect_equal(meanplot$labels$colour, "MS_batch")
  expect_equal(meanplot$labels$x, "order")
  expect_equal(meanplot$labels$y, "Mean_Intensity")
  
})


test_that("boxplot_plots", {
  data(example_proteome, package="proBatch")
  data(example_sample_annotation, package="proBatch")
  
  proteome <- example_proteome[1:20, ]
  expect_warning(boxplot <- plot_boxplot(proteome, example_sample_annotation, 
                                         batch_col = "MS_batch"))
  
  expect_equal(boxplot$labels$fill, "MS_batch")
  expect_equal(boxplot$label$group, "order")
  expect_equal(boxplot$label$x, "order")
  expect_equal(boxplot$label$y, "Intensity")
  
  expect_equal(boxplot$plot_env$color_by_batch, TRUE)
  expect_equal(boxplot$plot_env$facet_col, NULL)
  
})



