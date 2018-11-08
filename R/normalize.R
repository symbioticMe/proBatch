#' Data normalization and batch adjustment methods
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param batch_col column in `sample_annotation` that should be used for
#'   batch comparison
#' @param measure_col if `df_long` is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @name normalize
NULL
#> NULL


#' Log transformation of the data, ensuring that the row and column names
#' are retained
#'
#' @param data_matrix raw data matrix (features in rows and samples
#'   in columns)
#'
#' @return `data_matrix`-size matrix, with columns log2 transformed
#' @export
#'
#' @examples
log_transform <- function(data_matrix){
  data_matrix_log2 = log2(data_matrix + 1) 
  return(data_matrix_log2)
}

#' Quantile normalization of the data, ensuring that the row and column names
#' are retained
#'
#' @param data_matrix log transformed data matrix (features in rows and samples
#'   in columns)
#'
#' @return `data_matrix`-size matrix, with columns quantile-normalized
#' @export
#'
#' @examples
quantile_normalize <- function(data_matrix){
  q_norm_proteome = preprocessCore::normalize.quantiles(data_matrix)
  colnames(q_norm_proteome) = colnames(data_matrix)
  rownames(q_norm_proteome) = rownames(data_matrix)
  return(q_norm_proteome)
}


#' Median normalization of the data (global)
#'
#' @name normalize
#'
#' @return
#' @export
#'
#' @examples
normalize_medians_global <- function(df_long,
                                     sample_id_col = 'FullRunName',
                                    measure_col = 'Intensity'){
  df_normalized = df_long  %>%
    group_by_at(vars(one_of(sample_id_col))) %>%
    mutate(median_run = median(UQ(sym(measure_col)), na.rm = T)) %>%
    ungroup()
  df_normalized = df_normalized %>%
    mutate(median_global = median(UQ(sym(measure_col)), na.rm = T)) %>%
    mutate(diff = median_global - median_run) %>%
    mutate(Intensity_normalized = UQ(sym(measure_col))+diff)
  return(df_normalized)
}

#' Batch correction method allows correction of continuous sigal drift within batch and 
#' discrete difference across batches. 
#'
#' @name normalize
#' @param data_matrix raw data matrix (features in rows and samples
#'   in columns)
#' @param normalizeFunc global batch normalization method (`quantile` or `MedianCentering`)
#' @param log whether to log transform data matrix before normalization (`NULL`, `2` or `10`)
#'
#' @return `data_matrix`-size matrix, with columns normalized 
#' @export
#'
#' @examples
normalize <- function(data_matrix, normalizeFunc = "quantile", log = NULL){
  if(!is.null(log)){
    if(log == 2){
      data_matrix = log_transform(data_matrix)
    } else if(log == 10){
      data_matrix = log10(data_matrix + 1) 
    } else {
      stop("Only base 2 and base 10 logarithms are available.")
    }
  }
  
  if(normalizeFunc == "quantile"){
    normalized_matrix = quantile_normalize(data_matrix)
  } else if(normalizeFunc == "medianCentering"){
    df_long = matrix_to_long(matrix, feature_id_col = 'peptide_group_label',
                             measure_col = 'Intensity', sample_id_col = 'FullRunName')
    normalized_matrix = normalize_medians_global(df_long, sample_id_col = 'FullRunName',  measure_col = 'Intensity')
  } else {
    stop("Only quantile and median centering normalization methods are available")
  }
  
  return(normalized_matrix)
}
