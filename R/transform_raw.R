#' @inheritParams proBatch
#' @param log_base base of the logarithm for transformation
#' @param offset small positive number to prevent 0 conversion to \code{-Inf}
#' @name transform_raw_data

#' Log transformation of the data
#'
#'
#' @return `log_transform_df()` returns \code{df_long}-size data frame, with \code{measure_col} log transformed;
#' with old value in another column called "old_intensity" if "intensity" 
#' was the value of \code{measure_col};
#' `log_transform_matrix()` returns \code{data_matrix} of same dimensions as input
#' 
#' 
#' @examples 
#' log_transformed_df <- log_transform(example_proteome)
#' 
#' log_transformed_matrix <- log_transform_matrix(example_proteome_matrix, log_base = 10, offset = 1)
#' 
#' @export
#' @rdname transform_raw_data
#'
log_transform_df <- function(df_long, log_base = 2, offset = 1,
                             measure_col = 'Intensity'){
  if(!is.null(log_base)){
    df_long = df_long %>%
      mutate(!!(paste('old', measure_col, sep = '_')) := !!(sym(measure_col))) %>%
      mutate(!!(sym(measure_col)) := log(!!(sym(measure_col)) + offset, base = log_base))
  } else {
    warning("Log base is NULL, returning the original data frea")
  }
  return(df_long)
}

#' 
#' @export
#' @rdname transform_raw_data
#'
log_transform_matrix <- function(data_matrix, log_base = 2, offset = 1){
  if(!is.null(log_base)){
    if(log_base == 2){
      data_matrix_log = log2(data_matrix + offset)
    }else if(log_base ==10){
      data_matrix_log = log10(data_matrix + offset)
    }else {
      data_matrix_log = log(data_matrix + offset, base = log_base)
    }
  }
  return(data_matrix_log)
}