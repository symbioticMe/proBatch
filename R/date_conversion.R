#' Convert data/time to POSIXct
#'
#' convert date/time column of sample_annotation to POSIX format required to
#' keep number-like behavior
#'
#' @inheritParams proBatch
#' @param time_column name of the column(s) where run date & time are specified.
#'   These will be used to determine the run order
#' @param new_time_column name of the new column to which date&time will be
#'   converted to
#' @param dateTimeFormat POSIX format of the date and time. 
#'   See \code{\link{as.POSIXct}} from base R for details
#' @param tz for time zone, 'GMT' by default
#' @param locale for locale, 'en_US.UTF-8' by default
#'
#' @return sample annotation file with a new column \code{new_time_column} with
#'   POSIX-formatted date
#'
#' @family date
#' @examples 
#' date_to_posix <- dates_to_posix(example_sample_annotation, 
#' time_column = c('RunDate','RunTime'),
#' new_time_column = 'DateTime_new', 
#' dateTimeFormat = c("%b_%d", "%H:%M:%S"))
#'
#' @export
#'
dates_to_posix <- function(sample_annotation,
                           time_column = c('RunDate','RunTime'),
                           new_time_column = 'DateTime',
                           dateTimeFormat = c("%b_%d", "%H:%M:%S"), 
                           tz = 'GMT', locale = 'en_US.UTF-8'){
  old_locale <- Sys.getlocale("LC_TIME")
  Sys.setlocale("LC_TIME", locale)                          

  if (length(time_column) == 1){
    if(is.null(new_time_column)) new_time_column = time_column
    time_col = as.character(sample_annotation[[time_column]])
    sample_annotation[[new_time_column]] = as.POSIXct(time_col ,
                                                      format=dateTimeFormat, 
                                                      tz = tz)
  }
  else {
    sample_annotation = sample_annotation %>%
      mutate(dateTime = paste(!!!syms(time_column), sep=" ")) %>%
      mutate(dateTime = as.POSIXct(dateTime,
                                   format = paste(dateTimeFormat, 
                                                  collapse = ' '), 
                                                  tz = tz)) %>%
      rename(!!new_time_column := dateTime)
  }
  Sys.setlocale("LC_TIME", old_locale)
  return(sample_annotation)
}




#' Convert date/time to POSIXct and rank samples by it
#'
#' Converts date/time columns fo sample_annotation to POSIXct format and
#' calculates sample run rank in order column
#'
#' @inheritParams dates_to_posix
#'
#' @param new_order_col name of column with generated the order of sample run 
#'  based on time columns
#' @param instrument_col column, denoting different 
#' instrument used for measurements
#' @return sample annotation file with a new column \code{new_time_column} with
#'   POSIX-formatted date & \code{new_order_col} used 
#'   in some diagnostic plots (e.g.
#'   \code{\link{plot_iRT}}, \code{\link{plot_sample_mean}})
#' @examples 
#' sample_annotation_wOrder <- date_to_sample_order(
#' example_sample_annotation,
#' time_column = c('RunDate','RunTime'),
#' new_time_column = 'new_DateTime',
#' dateTimeFormat = c("%b_%d", "%H:%M:%S"),
#' new_order_col = 'new_order',
#' instrument_col = NULL)
#' 
#' @export
#'
#' @name date_to_sample_order
date_to_sample_order <- function(sample_annotation,
                                 time_column = c('RunDate','RunTime'),
                                 new_time_column = 'DateTime',
                                 dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                                 new_order_col = 'order',
                                 instrument_col = 'instrument'){
  sample_annotation = dates_to_posix(sample_annotation = sample_annotation,
                                     time_column = time_column,
                                     new_time_column = new_time_column,
                                     dateTimeFormat = dateTimeFormat)
  sample_annotation = sample_annotation %>% arrange(!!(sym(new_time_column)))
  if (!is.null(instrument_col)){
    sample_annotation = sample_annotation %>%
      group_by_at(vars(one_of(instrument_col))) %>%
      mutate(!!(sym(new_order_col)) := rank(!!sym(new_time_column))) %>%
      ungroup()
  } else {
    sample_annotation[[new_order_col]] = 
      rank(sample_annotation[[new_time_column]])
  }
  return(sample_annotation)
}
