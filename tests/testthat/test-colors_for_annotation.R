context("colors_for_annotation")


test_that("map_factors_to_colors", {
  data(example_sample_annotation, package="proBatch")

  factor_columns = c('MS_batch', "Strain")
  sample_annt <- example_sample_annotation[,factor_columns]
  
  expect_warning(sample_color_annt <- map_factors_to_colors(sample_annt))
  
  expect_equal(names(sample_color_annt$MS_batch[1]), "Batch_1")
  expect_equal(sample_color_annt$MS_batch[[1]], "turquoise")
  expect_equal(sample_color_annt$MS_batch[[2]], "blue")
})



test_that("map_numbers_to_colors", {
  data(example_sample_annotation, package="proBatch")
  
  numeric_columns = c("DateTime", "order")
  sample_annt <- example_sample_annotation[,numeric_columns]
  sample_color_annt <- map_numbers_to_colors(sample_annt)
  
  expect_equal(names(sample_color_annt$color_list), c("DateTime", "order"))
  expect_equal(sample_color_annt$color_list$DateTime[[1]],"#A6611A")
  expect_equal(sample_color_annt$color_list$order[[1]], "#D01C8B")
  
  expect_equal(names(sample_color_annt$color_list$order[1]), "(0.768,24.2]")
  expect_equal(names(sample_color_annt$color_list$order[2]), "(24.2,47.4]")

})




test_that("sample_annotation_to_colors", {
  data(example_sample_annotation, package="proBatch")
  
  expect_warning(color_scheme <- sample_annotation_to_colors(example_sample_annotation, 
                                               factor_columns = c('MS_batch','EarTag', "Strain", 
                                                                  "Diet", "digestion_batch", "Sex"),
                                               date_columns = 'DateTime',
                                               numeric_columns = c('order')))
  
  expect_equal(names(color_scheme), c("list_of_colors", "color_df", "sample_annotation"))
  expect_equal(names(color_scheme$list_of_colors), c("MS_batch", "EarTag", "Strain", "Diet",  "digestion_batch", "Sex",            
                                        "FullRunName", "RunDate", "RunTime", "DateTime", "order"))
  
  expect_equal(color_scheme$list_of_colors$MS_batch[[1]], "turquoise")
  expect_equal(color_scheme$list_of_colors$MS_batch[[2]],  "blue")
  
  expect_equal(names(color_scheme$list_of_colors$Sex)[1], "F")
  expect_equal(names(color_scheme$list_of_colors$Sex)[2], "M")
  
})


