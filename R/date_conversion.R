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
dates_to_posix <- function(sample_annotation,
                           time_column = c('RunDate','RunTime'),
                           new_time_column = NULL,
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
#' @import dplyr
#' @import rlang
#' @export
#'
#' @examples
date_to_sample_order <- function(sample_annotation,
                                 time_column = c('RunDate','RunTime'),
                                 new_time_column = 'DateTime',
                                 dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                 order_column = 'order',
                                 instrument_col = 'instrument'){
  sample_annotation = dates_to_posix(sample_annotation = sample_annotation,
                                     time_column = time_column,
                                     new_time_column = new_time_column,
                                     dateTimeFormat = dateTimeFormat)
  sample_annotation = sample_annotation %>% arrange(UQ(rlang::sym(new_time_column)))
  if (!is.null(instrument_col)){
    sample_annotation = sample_annotation %>%
      group_by_at(vars(one_of(instrument_col))) %>%
      mutate(UQ(rlang::sym(order_column)) := rank(!!rlang::sym(new_time_column))) %>%
      ungroup()
  } else {
    sample_annotation[[order_column]] = rank(sample_annotation[[new_time_column]])
  }
  return(sample_annotation)
}

#' Identify stretches of time between runs that are long and split a batches by them
#'
#' @param date_vector POSIX or numeric-like vector corresponding to the sample
#' MS profile acquisition timepoint
#' @param threshold time difference that would mean there was an interruption
#' @param minimal_batch_size minimal number of samples in a batch
#' @param batch_name string with a self-explanatory name for the batch (e.g. `MS_batch` for MS-proteomics) to which batch number will be added
#'
#' @return vector of batches for each sample
#' @export
#'
#' @examples
define_batches_by_MS_pauses_within_instrument <- function(date_vector, threshold,
                                        minimal_batch_size = 5,
                                        batch_name = 'MS_batch'){
  diff = diff(date_vector)
  tipping_points = which(diff > threshold)
  batch_size = diff(tipping_points)
  if (any(batch_size <= minimal_batch_size)){
    warning('some batches are too small, merging with the previous')
    batch_correction = ifelse(batch_size <= minimal_batch_size, batch_size, 0)
    tipping_points = unique(tipping_points+ c(batch_correction, 0))
  }
  tipping_points = c(0, tipping_points, length(date_vector))
  batch_idx = rep(1:(length(tipping_points) -1),
                  times = diff(tipping_points))
  batch_ids = paste(batch_name, batch_idx, sep = '_')
  return(batch_ids)
}

#' convert by batch
#'
#' @param sample_annotation
#' @param threshold
#' @param runtime_col
#' @param minimal_batch_size
#' @param instrument_col
#' @param batch_name
#'
#' @return
#' @export
#'
#' @examples
define_batches_by_MS_pauses <- function(sample_annotation,
                                        threshold,
                                        runtime_col = 'RunDateTime',
                                        minimal_batch_size = 5,
                                        instrument_col = 'instr',
                                        batch_name = 'MS_batch'){

  if (!is.null(instrument_col)){
    #batch_name = paste(sample_annotation[[instrument_col]], batch_name, sep = ':')
    sample_annotation = sample_annotation %>%
      mutate(batch_name = paste(!!rlang::sym(instrument_col), batch_name, sep = ':')) %>%
      arrange(UQ(rlang::sym(instrument_col)), UQ(rlang::sym(runtime_col))) %>%
      group_by_at(vars(one_of(instrument_col))) %>%
      mutate(batch_id = define_batches_by_MS_pauses_within_instrument(!!rlang::sym(runtime_col),
                                                                        threshold = threshold,
                                                                        minimal_batch_size = minimal_batch_size,
                                                                        batch_name = batch_name))%>%
        rename(!!batch_name := batch_id)

  } else {
    sample_annotation = sample_annotation
      mutate(batch_id = define_batches_by_MS_pauses_within_instrument(!!rlang::sym(runtime_col),
                                                                        threshold = threshold,
                                                                        minimal_batch_size = minimal_batch_size,
                                                                        batch_name = batch_name))%>%
        rename(!!batch_name := batch_id)
  }
  return(sample_annotation)
}
