#' Convert data/time to POSIXct
#'
#' convert date/time column of sample_annotation to POSIX format required to
#' keep number-like behaviour
#'
#' @inheritParams proBatch
#' @param time_column name of the column(s) where run date & time are specified.
#'   These will be used to determine the run order
#' @param new_time_column name of the new column to which date&time will be
#'   converted to
#' @param dateTimeFormat POSIX format of the date and time. See \code{\link{as.POSIXct}}
#'   from base R for details
#'
#' @return sample annotation file with column names as 'new_time_column' with
#'   POSIX-formatted date
#'
#' @family date
#'
#' @export
#'
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
      mutate(dateTime = paste(!!!syms(time_column), sep=" ")) %>%
      mutate(dateTime = lubridate::as.POSIXct(dateTime,
                                   format = paste(dateTimeFormat, collapse = ' '))) %>%
      rename(!!new_time_column := dateTime)
  }
  return(sample_annotation)
}


#' Convert date/time to POSIXct and rank samples by it
#'
#' Converts date/time columns fo sample_annotation to POSIXct format and
#' calculates sample run rank in order column
#'
#' @inheritParams dates_to_posix
#'
#' @return sample annotation file with column names as 'new_time_column' with
#'   POSIX-formatted date & \code{order_col} used in some diagnostic plots (e.g.
#'   \code{\link{plot_iRTs}}, \code{\link{plot_sample_mean}})
#'
#' @export
#'
date_to_sample_order <- function(sample_annotation,
                                 time_column = c('RunDate','RunTime'),
                                 new_time_column = 'DateTime',
                                 dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                 order_col = 'order',
                                 instrument_col = 'instrument'){
  sample_annotation = dates_to_posix(sample_annotation = sample_annotation,
                                     time_column = time_column,
                                     new_time_column = new_time_column,
                                     dateTimeFormat = dateTimeFormat)
  sample_annotation = sample_annotation %>% arrange(UQ(sym(new_time_column)))
  if (!is.null(instrument_col)){
    sample_annotation = sample_annotation %>%
      group_by_at(vars(one_of(instrument_col))) %>%
      mutate(UQ(sym(order_col)) := rank(!!sym(new_time_column))) %>%
      ungroup()
  } else {
    sample_annotation[[order_col]] = rank(sample_annotation[[new_time_column]])
  }
  return(sample_annotation)
}


#' Batch by date/time
#'
#' Identify long stretches of time between samples and split them into batches.
#' Most users are going to want to call define_batches_by_MS_pauses, rather then
#' this function
#'
#' @param date_vector POSIX or numeric-like vector corresponding to the sample
#'   MS profile acquisition timepoint
#' @param threshold time difference that would mean there was an interruption
#' @param minimal_batch_size minimal number of samples in a batch
#' @param batch_name string with a self-explanatory name for the batch (e.g.
#'   `MS_batch` for MS-proteomics) to which batch number will be added
#'
#' @return vector of batches for each sample
#'
#' @export
#'
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


#' Batch by date/time and instrument
#'
#' Identify long stretches of time between samples on a per instrument basis and
#' split them into batches.
#'
#'
#' @inheritParams proBatch
#' @inheritParams define_batches_by_MS_pauses_within_instrument
#' @param runtime_col POSIX or numeric-like column corresponding to the sample
#'   MS profile acquisition timepoint
#' @param instrument_col column specifying MS instrument used to acquired data
#'   (to account for the presence of multiple instruments in
#'   \code{sample_annotation})
#'
#' @return \code{sample_annotation} data matrix with an additional
#'   column to indicate sample batching by MS run time and instrument.
#'
#' @export
#'
define_batches_by_MS_pauses <- function(sample_annotation,
                                        threshold,
                                        runtime_col = 'RunDateTime',
                                        minimal_batch_size = 5,
                                        instrument_col = 'instr',
                                        batch_name = 'MS_batch'){

  if (!is.null(instrument_col)){
    sample_annotation = sample_annotation %>%
      mutate(batch_name = paste(!!sym(instrument_col), batch_name, sep = ':')) %>%
      arrange(UQ(sym(instrument_col)), UQ(sym(runtime_col))) %>%
      group_by_at(vars(one_of(instrument_col))) %>%
      mutate(batch_id = define_batches_by_MS_pauses_within_instrument(!!sym(runtime_col),
                                                                        threshold = threshold,
                                                                        minimal_batch_size = minimal_batch_size,
                                                                        batch_name = batch_name))%>%
        rename(!!batch_name := batch_id)

  } else {
    sample_annotation = sample_annotation
      mutate(batch_id = define_batches_by_MS_pauses_within_instrument(!!sym(runtime_col),
                                                                        threshold = threshold,
                                                                        minimal_batch_size = minimal_batch_size,
                                                                        batch_name = batch_name))%>%
        rename(!!batch_name := batch_id)
  }
  return(sample_annotation)
}
