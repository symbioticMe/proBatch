#' Data normalization methods
#' 
#' @description Normalization of raw (usually log-transformed) data. 
#' Normalization brings the samples to the same scale.
#' Currently the following normalization functions are implemented:
#' #' \enumerate{
#'   \item Quantile normalization: `quantile_normalize_dm()`. 
#'   Quantile normalization of the data.
#'   \item Median normalization: `normalize_sample_medians_dm()`. 
#'   Normalization by centering sample medians to global median of the data
#' }
#' Alternatively, one can call normalization function with `normalize_data_dm()` 
#' wrapper.
#' 
#' 
#' @inheritParams proBatch
#' @param normalize_func global batch normalization method 
#' (`quantile` or `MedianCentering`)
#' @param log_base whether to log transform data matrix 
#' before normalization (e.g. `NULL`, `2` or `10`)
#' @param offset small positive number to prevent 0 conversion to \code{-Inf}
#' 
#' @return the data in the same format as input (\code{data_matrix} or 
#' \code{df_long}).
#' For \code{df_long} the data frame stores the original values of 
#' \code{measure_col}
#' in another column called "preNorm_intensity" if "intensity", and the 
#' normalized values in \code{measure_col} column.
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
                                  measure_col = 'Intensity',
                                  no_fit_imputed = TRUE,
                                  qual_col = NULL,
                                  qual_value = 2,
                                  keep_all = 'default'){
  
  if(is.null(qual_col) & no_fit_imputed){
    warning('imputed value flag column is NULL, changing no_fit_imputed to FALSE')
    no_fit_imputed = FALSE
  }
  
  if(no_fit_imputed){
    if(!(qual_col %in% names(df_long))){
      stop("imputed value flag column (qual_col) is not in the data frame!")
    }
    message('removing imputed values (requants) from the matrix')
    data_matrix = long_to_matrix(df_long, 
                                 feature_id_col = feature_id_col, 
                                 measure_col = measure_col,
                                 sample_id_col = sample_id_col, 
                                 qual_col = qual_col, 
                                 qual_value = qual_value)
  } else {
    if(!is.null(qual_col) && (qual_col %in% names(df_long))){
      warning('imputed value (requant) column is in the data, are you sure you
                want to use imputed (requant) values in quantile inference?')
    }
    data_matrix = long_to_matrix(df_long, 
                                 feature_id_col = feature_id_col, 
                                 measure_col = measure_col,
                                 sample_id_col = sample_id_col, 
                                 qual_col = NULL)
  }
  
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
    merge(df_long %>% select(-one_of(setdiff(names(normalized_df), 
                                             c(feature_id_col, sample_id_col, measure_col)))), 
          by = c(feature_id_col, sample_id_col))
  
  if(!is.null(qual_col) && qual_col %in% names(normalized_df)){
    normalized_df = switch (keep_all,
                           all = normalized_df,
                           default = normalized_df,
                           minimal = normalized_df %>%
                             dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                             old_measure_col, qual_col, qual_value)))
    )
  } else {
    normalized_df = switch (keep_all,
                            all = normalized_df,
                            default = normalized_df, 
                            minimal = normalized_df %>%
                              dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                              old_measure_col))))
  }
  
  return(normalized_df)
}

#' 
#' @export
#' @rdname normalize
#'
normalize_sample_medians_dm <- function(data_matrix){
  df_long = matrix_to_long(data_matrix, 
                           feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity', 
                           sample_id_col = 'FullRunName')
  normalized_df = normalize_sample_medians_df(df_long, 
                                              sample_id_col = 'FullRunName', 
                                              measure_col = 'Intensity')
  normalized_matrix = long_to_matrix(normalized_df)
  return(normalized_matrix)
}

#' 
#' @export
#' @rdname normalize
#'
normalize_sample_medians_df <- function(df_long,
                                        feature_id_col = 'peptide_group_label',
                                        sample_id_col = 'FullRunName',
                                        measure_col = 'Intensity',
                                        no_fit_imputed = FALSE,
                                        qual_col = NULL,
                                        qual_value = 2,
                                        keep_all = 'default'){
  if(no_fit_imputed){
    if(!(qual_col %in% names(df_long))){
      stop("imputed value flag column (qual_col) is not in the data frame!")
    }
    message('removing imputed values (requants) from the matrix')
    df_long = df_long %>%
      mutate(!!sym(measure_col) := ifelse(!!sym(qual_col) == qual_value, 
                                          NA, measure_col))
  } else {
    if(!is.null(qual_col) && (qual_col %in% names(df_long))){
      warning('imputed value (requant) column is in the data, are you sure you
              want to use imputed (requant) values in sample median inference?')
    }
  }
  
  normalized_df = df_long  %>%
    group_by_at(vars(one_of(sample_id_col))) %>%
    mutate(median_run = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup()
  
  old_measure_col = paste('preNorm', measure_col, sep = '_')
  
  normalized_df = normalized_df %>%
    mutate(median_global = median(!!(sym(measure_col)), na.rm = TRUE),
           
           !!(old_measure_col) := !!(sym(measure_col))) %>%
    mutate(diff_norm = median_global - median_run) %>%
    mutate(!!(sym(measure_col)) := !!(sym(measure_col))+diff_norm)
  
  if(!is.null(qual_col) && (qual_col %in% names(normalized_df))){
    normalized_df = switch (keep_all,
                            all = normalized_df,
                            default = normalized_df,
                            minimal = normalized_df %>%
                              dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                              old_measure_col, qual_col, qual_value)))
    )
  } else {
    normalized_df = switch (keep_all,
                            all = normalized_df,
                            default = normalized_df, 
                            minimal = normalized_df %>%
                              dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                              old_measure_col))))
  }
  
  return(normalized_df)
}

#' 
#' @export
#' @rdname normalize
#' 
#' 
normalize_data_dm <- function(data_matrix, 
                              normalize_func = c("quantile", "medianCentering"), 
                           log_base = NULL, offset = 1){
  
  normalize_func <- match.arg(normalize_func)    
  if(!is.null(log_base)){
    data_matrix = log_transform_dm(data_matrix, log_base = log_base)
  }
  
  if(normalize_func == "quantile"){
    normalized_matrix = quantile_normalize_dm(data_matrix)
  } else if(normalize_func == "medianCentering"){
    normalized_matrix = normalize_sample_medians_dm(data_matrix)
  } else {
    stop("Only quantile and median centering normalization methods implemented")
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
                              measure_col = 'Intensity',
                              no_fit_imputed = TRUE,
                              qual_col = NULL,
                              qual_value = 2,
                              keep_all = 'default'){
  
  normalize_func <- match.arg(normalize_func)
  
  
  if(!is.null(log_base)){
    df_long = log_transform_df(df_long, log_base = log_base, offset = offset)
  }
  
  if(normalize_func == "quantile"){
    normalized_df = quantile_normalize_df(df_long, 
                                          feature_id_col = feature_id_col, 
                                          measure_col = measure_col,
                                          sample_id_col = sample_id_col,
                                          keep_all = keep_all, 
                                          no_fit_imputed = no_fit_imputed,
                                          qual_col = qual_col,
                                          qual_value = qual_value)
  } else if(normalize_func == "medianCentering"){
    normalized_df = normalize_sample_medians_df(df_long, 
                                                sample_id_col = sample_id_col, 
                                                measure_col = measure_col,
                                                keep_all = keep_all, 
                                                no_fit_imputed = no_fit_imputed,
                                                qual_col = qual_col,
                                                qual_value = qual_value)
  } else {
    stop("Only quantile and median centering normalization methods implemented")
  }
  
  return(normalized_df)
}
