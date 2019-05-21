#' Data normalization methods
#' 
#' @description Normalization of raw (usually log-transformed) data. Normalization brings the samples to the same scale
#' Currently the following normalization functions are implemented:
#' #' \enumerate{
#'   \item Quantile normalization: `quantile_normalize()`. Quantile normalization of the data.
#'   \item Median normalization: `normalize_sample_medians()`. Normalization by centering sample medians to global median of the data
#' }
#' Alternatively, one can call normalization function with `normalize_data()` wrapper.
#' 
#' 
#' @inheritParams proBatch
#' @param normalize_func global batch normalization method 
#' (`quantile` or `MedianCentering`)
#' @param log_base whether to log transform data matrix 
#' before normalization (`NULL`, `2` or `10`)
#' 
#' @return the data in the same format as input (\code{data_matrix} or \code{df_long}).
#' For \code{df_long} the data frame stores the original values of \code{measure_col}
#' in another column called "pre_norm_intensity" if "intensity", and the normalized values
#' in \code{measure_col} column.
#' 
#' @examples 
#' 
#' #Quantile normalization:
#' quantile_normalized_matrix <- quantile_normalize(example_proteome_matrix)
#' 
#' #Median centering:
#' median_normalized_df <- normalize_sample_medians(example_proteome)
#' 
#' #Transform the data in one go:
#' quantile_normalized_matrix <- normalize_data(example_proteome_matrix, 
#' normalizeFunc = "quantile", log_base = 2)
#' 
#' @name normalize
NULL

#' 
#' @export
#' @rdname normalize
#'
quantile_normalize <- function(data_matrix){
  q_norm_proteome = normalize.quantiles(data_matrix)
  colnames(q_norm_proteome) = colnames(data_matrix)
  rownames(q_norm_proteome) = rownames(data_matrix)
  return(q_norm_proteome)
}

#' 
#' @export
#' @rdname normalize
#'
normalize_sample_medians <- function(df_long,
                                     sample_id_col = 'FullRunName',
                                     measure_col = 'Intensity'){
  df_normalized = df_long  %>%
    group_by_at(vars(one_of(sample_id_col))) %>%
    mutate(median_run = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup()
  df_normalized = df_normalized %>%
    mutate(median_global = median(!!(sym(measure_col)), na.rm = TRUE),
           !!(paste('pre_norm', measure_col, sep = '_')) := !!(sym(measure_col))) %>%
    mutate(diff = median_global - median_run) %>%
    mutate(!!(sym(measure_col)) := !!(sym(measure_col))+diff)
  return(df_normalized)
}

#' 
#' @export
#' @rdname normalize
#' 
#' 
normalize_data <- function(data_matrix, normalize_func = c("quantile", "medianCentering"), 
                           log_base = NULL){
  
  normalize_func <- match.arg(normalize_func)    
  if(!is.null(log_base)){
    data_matrix = log_transform_matrix(data_matrix, log_base = log_base)
  }
  
  if(normalize_func == "quantile"){
    normalized_matrix = quantile_normalize(data_matrix)
  } else if(normalize_func == "medianCentering"){
    df_long = matrix_to_long(data_matrix, 
                             feature_id_col = 'peptide_group_label',
                             measure_col = 'Intensity', 
                             sample_id_col = 'FullRunName')
    normalized_df = normalize_sample_medians(df_long, 
                                             sample_id_col = 'FullRunName', 
                                             measure_col = 'Intensity')
    normalized_matrix = long_to_matrix(normalized_df)
  } else {
    stop("Only quantile and median centering normalization methods are available")
  }
  
  return(normalized_matrix)
}
