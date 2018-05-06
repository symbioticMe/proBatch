#' convert date/time column of sample_annotation to POSIX format
#' required to keep number-like behaviour
#'
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names)
#' 2) biological and 3) technical covariates (batches etc)
#' @param time_column name of the column(s) where run date & time are specified.
#' These will be used to determine the run order
#' @param new_time_column name of the new column to which date&time will be converted to
#' @param dateTimeFormat POSIX format of the date and time. See `as.POSIXct` from base R for details
#'
#' @return sample annotation file with column names as 'new_time_column' with POSIX-formatted date
#' @export
#'
#' @examples
dates_to_posix <- function(sample_annotation, time_column, new_time_column = NULL,
                           dateTimeFormat = c("%b_%d", "%H:%M:%S")){
  if (length(time_column) == 1){
    if(is.null(new_time_column)) new_time_column = time_column
    time_col = as.character(sample_annotation[[time_column]])
    sample_annotation[[new_time_column]] = as.POSIXct(time_col ,
                                                      format=dateTimeFormat)
  }
  else {
    sample_annotation = sample_annotation %>%
      mutate(dateTime = paste(!!!rlang::syms(time_column), sep=" ")) %>%
      mutate(dateTime = as.POSIXct(dateTime,
                                   format = paste(dateTimeFormat, collapse = ' '))) %>%
      rename(!!new_time_column := dateTime)
  }
  return(sample_annotation)
}

#' convert date to order
#'
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names)
#' 2) biological and 3) technical covariates (batches etc)
#' @param time_column name of the column(s) where run date & time are specified.
#' These will be used to determine the run order
#' @param new_time_column name of the new column to which date&time will be converted to
#' @param dateTimeFormat POSIX format of the date and time. See `as.POSIXct` from base R for details
#' @param order_column name of the new column that determines sample order.
#' Will be used for certain diagnostics and normalisations
#'
#' @return sample annotation file with column names as 'new_time_column'
#' with POSIX-formatted date & `order_column` used in some diagnostic plots (`plot_iRTs`, `plot_sample_mean`)
#' @export
#'
#' @examples
date_to_sample_order <- function(sample_annotation, time_column,
                                 new_time_column = 'DateTime',
                                 dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                 order_column = 'order'){
  sample_annotation = dates_to_posix(sample_annotation = sample_annotation,
                                     time_column = time_column,
                                     new_time_column = new_time_column,
                                     dateTimeFormat = dateTimeFormat)
  sample_annotation[[order_column]] = rank(sample_annotation[[new_time_column]])
  return(sample_annotation)
}
