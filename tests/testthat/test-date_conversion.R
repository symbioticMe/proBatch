context("date_conversion")


test_that("dates_to_posix", {
  data(example_sample_annotation, package="proBatch")

  sample_annotation <- example_sample_annotation[1:5,]
  sample_annotation$DateTime <- NULL
  
  new_annotation <- dates_to_posix(sample_annotation, 
                                   time_column = c('RunDate','RunTime'),
                                   new_time_column = 'DateTime', 
                                   dateTimeFormat = c("%b_%d", "%H:%M:%S"))
  
  expect_equal(new_annotation$RunDate[1], "Oct_05")
  expect_equal(as.character(new_annotation$RunTime[1]), "18:35:00")
  expect_equal(new_annotation$DateTime[1], as.POSIXlt("2019-10-05 18:35:00 CEST"))
  
  expect_is(new_annotation$DateTime, "POSIXct")
  expect_is(new_annotation$RunTime, "difftime")
})


test_that("date_to_sample_order", {
  data(example_sample_annotation, package="proBatch")
  
  sample_annotation_test <- example_sample_annotation[1:5,]
  sample_annotation_test$DateTime <- NULL
  sample_annotation_test$order <- NULL
  
  new_annotation_worder <- date_to_sample_order(sample_annotation_test,
                                          time_column = c('RunDate','RunTime'),
                                          new_time_column = 'new_DateTime',
                                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                          new_order_col = 'new_order',
                                          instrument_col = NULL)
  
  expect_equal(new_annotation_worder$RunDate[1], "Oct_05")
  expect_equal(as.character(new_annotation_worder$RunTime[1]), "18:35:00")
  
  expect_is(new_annotation_worder$new_DateTime, "POSIXct")
  expect_is(new_annotation_worder$RunTime, "difftime")
  expect_length(new_annotation_worder$new_order, nrow(sample_annotation_test))
  
  expect_equal(new_annotation_worder$new_order[1], 1)
  expect_equal(new_annotation_worder$new_order[2], 2)
  
  

})
