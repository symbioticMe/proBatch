#' Data normalization methods
#' 
#' @description Normalization of raw (usually log-transformed) data. 
#' Normalization brings the samples to the same scale.
#' Currently the following normalization functions are implemented:
#' #' \enumerate{
#'   \item Quantile normalization: `quantile_normalize()`. 
#'   Quantile normalization of the data.
#'   \item Median normalization: `normalize_sample_medians()`. 
#'   Normalization by centering sample medians to global median of the data
#' }
#' Alternatively, one can call normalization function with `normalize_data()` 
#' wrapper.
#' 
#' 
#' @inheritParams proBatch
#' @param normalize_func global batch normalization method 
#' (`quantile` or `MedianCentering`)
#' @param log_base whether to log transform data matrix 
#' before normalization (e.g. `NULL`, `2` or `10`)
#' 
#' @return the data in the same format as input (\code{data_matrix} or 
#' \code{df_long}).
#' For \code{df_long} the data frame stores the original values of 
#' \code{measure_col}
#' in another column called "preNorm_intensity" if "intensity", and the 
#' normalized values
#' in \code{measure_col} column.
#' 
#' @examples 
#' 
#' #Quantile normalization:
#' quantile_normalized_matrix <- quantile_normalize_dm(example_proteome_matrix)
#' 
#' #Median centering:
#' median_normalized_df <- normalize_sample_medians_df(example_proteome)
#' 
#' #Transform the data in one go:
#' quantile_normalized_matrix <- normalize_data_dm(example_proteome_matrix, 
#' normalize_func = "quantile", log_base = 2, offset = 1)
#' 
#' @name normalize
NULL

#' 
#' @export
#' @rdname normalize
#'
quantile_normalize_dm <- function(data_matrix){
  q_norm_proteome = normalize.quantiles(data_matrix)
  colnames(q_norm_proteome) = colnames(data_matrix)
  rownames(q_norm_proteome) = rownames(data_matrix)
  return(q_norm_proteome)
}

#' 
#' @export
#' @rdname normalize
#'
quantile_normalize_df <- function(df_long,
                                  feature_id_col = 'peptide_group_label',
                                  sample_id_col = 'FullRunName', 
                                  measure_col = 'Intensity'){
  
  
  data_matrix = long_to_matrix(df_long, 
                               feature_id_col = feature_id_col, 
                               measure_col = measure_col,
                               sample_id_col = sample_id_col)
  q_norm_proteome = normalize.quantiles(data_matrix)
  colnames(q_norm_proteome) = colnames(data_matrix)
  rownames(q_norm_proteome) = rownames(data_matrix)
  normalized_df = matrix_to_long(q_norm_proteome,
                                 feature_id_col = feature_id_col, 
                                 measure_col = measure_col,
                                 sample_id_col = sample_id_col)
  
  old_measure_col = paste('preNorm', measure_col, sep = '_')
  
  df_long = df_long %>% 
    rename(!!(old_measure_col) := !!(sym(measure_col))) 
  
  normalized_df = normalized_df %>%
    merge(df_long, by = c(feature_id_col, sample_id_col))
  return(normalized_df)
}

#' 
#' @export
#' @rdname normalize
#'
normalize_sample_medians_dm <- function(data_matrix,
                                        sample_id_col = 'FullRunName',
                                        measure_col = 'Intensity'){
  df_normalized = df_long  %>%
    group_by_at(vars(one_of(sample_id_col))) %>%
    mutate(median_run = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup()
  
  old_measure_col = paste('preNorm', measure_col, sep = '_')
  
  df_normalized = df_normalized %>%
    mutate(median_global = median(!!(sym(measure_col)), na.rm = TRUE),
           
           !!(old_measure_col) := !!(sym(measure_col))) %>%
    mutate(diff = median_global - median_run) %>%
    mutate(!!(sym(measure_col)) := !!(sym(measure_col))+diff)
  return(df_normalized)
}

#' 
#' @export
#' @rdname normalize
#'
normalize_sample_medians_df <- function(df_long,
                                     sample_id_col = 'FullRunName',
                                     measure_col = 'Intensity'){
  df_normalized = df_long  %>%
    group_by_at(vars(one_of(sample_id_col))) %>%
    mutate(median_run = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup()
  
  old_measure_col = paste('preNorm', measure_col, sep = '_')
  
  df_normalized = df_normalized %>%
    mutate(median_global = median(!!(sym(measure_col)), na.rm = TRUE),
           
           !!(old_measure_col) := !!(sym(measure_col))) %>%
    mutate(diff = median_global - median_run) %>%
    mutate(!!(sym(measure_col)) := !!(sym(measure_col))+diff)
  return(df_normalized)
}

#' 
#' @export
#' @rdname normalize
#' 
#' 
normalize_data_dm <- function(data_matrix, normalize_func = c("quantile", "medianCentering"), 
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

#' @export
#' @rdname normalize
#' 
#' 
normalize_data_df <- function(df_long, 
                              normalize_func = c("quantile", "medianCentering"), 
                              log_base = NULL, offset = 1,
                              feature_id_col = 'peptide_group_label',
                              sample_id_col = 'FullRunName', 
                              measure_col = 'Intensity'){
  
  normalize_func <- match.arg(normalize_func)
  
  
  if(!is.null(log_base)){
    data_matrix = log_transform_df(df_long, log_base = log_base, offset = offset)
  }
  
  if(normalize_func == "quantile"){
    data_matrix = long_to_matrix(df_long, 
                                 feature_id_col = feature_id_col, 
                                 measure_col = measure_col,
                                 sample_id_col = sample_id_col)
    normalized_matrix = quantile_normalize(data_matrix)
    normalized_df = matrix_to_long(normalized_matrix,
                                   feature_id_col = feature_id_col, 
                                   measure_col = measure_col,
                                   sample_id_col = sample_id_col)
  } else if(normalize_func == "medianCentering"){
    normalized_df = normalize_sample_medians_df(df_long, 
                                             sample_id_col = 'FullRunName', 
                                             measure_col = 'Intensity')
  } else {
    stop("Only quantile and median centering normalization methods are available")
  }
  
  return(normalized_df)
}
