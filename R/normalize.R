#' Data normalization methods
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if `df_long` is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @name normalize
NULL


#' Log transformation of the data, ensuring that the row and column names
#' are retained
#'
#' @param data_matrix raw data matrix (features in rows and samples
#'   in columns)
#' @param log_base base of the logarithm for transformation
#'
#' @return `data_matrix`-size matrix, with columns log2 transformed
#' @examples 
#' log_transformed_matrix <- log_transform(example_proteome_matrix)
#' 
#' @export
#'
log_transform <- function(data_matrix, log_base = 2){
  if(!is.null(log_base)){
    if(log_base == 2){
      data_matrix_log = log2(data_matrix + 1)
    }else if(log_base ==10){
      data_matrix_log = log10(data_matrix + 1)
    }else {
      stop("Only log_base = 2 or log_base = 10 is applicable")
    }
  }
    return(data_matrix_log)
}

#' Quantile normalization of the data, ensuring that the row and column names
#' are retained
#'
#' @param data_matrix log transformed data matrix (features in rows and samples
#'   in columns)
#'
#' @return `data_matrix`-size matrix, with columns quantile-normalized
#' @examples 
#' quantile_normalized_matrix = quantile_normalize(example_proteome_matrix)
#' 
#' @export
#'
quantile_normalize <- function(data_matrix){
    q_norm_proteome = preprocessCore::normalize.quantiles(data_matrix)
    colnames(q_norm_proteome) = colnames(data_matrix)
    rownames(q_norm_proteome) = rownames(data_matrix)
    return(q_norm_proteome)
}


#' Normalization by centering sample medians to global median of the data
#'
#' @param df_long log transformed long format data matrix (see `df_long`)
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if `df_long` is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @return `df_long`-size matrix, with intensity scaled to global median
#' @examples 
#' median_normalized_matrix = normalize_sample_medians(example_proteome)
#' 
#' @export
#'
normalize_sample_medians <- function(df_long,
                                     sample_id_col = 'FullRunName',
                                     measure_col = 'Intensity'){
    df_normalized = df_long  %>%
        group_by_at(vars(one_of(sample_id_col))) %>%
        mutate(median_run = median(UQ(sym(measure_col)), na.rm = TRUE)) %>%
        ungroup()
    df_normalized = df_normalized %>%
        mutate(median_global = median(UQ(sym(measure_col)), na.rm = TRUE)) %>%
        mutate(diff = median_global - median_run) %>%
        mutate(Intensity_normalized = UQ(sym(measure_col))+diff)
    return(df_normalized)
}

#' Normalization brings the samples to the same scale
#'
#' @param data_matrix raw data matrix (features in rows and samples
#'   in columns)
#' @param normalizeFunc global batch normalization method 
#' (`quantile` or `MedianCentering`)
#' @param log_base whether to log transform data matrix 
#' before normalization (`NULL`, `2` or `10`)
#'
#' @return `data_matrix`-size matrix, with columns normalized 
#' @export
#'
#' @examples
#' quantile_normalized_matrix = normalize_data(example_proteome_matrix, 
#' normalizeFunc = "quantile", log_base = 2)
#' 
normalize_data <- function(data_matrix, normalizeFunc = c("quantile", "medianCentering"), 
                           log_base = NULL){
  
    normalizeFunc <- match.arg(normalizeFunc)    
    if(!is.null(log_base)){
        data_matrix = log_transform(data_matrix, log_base = log_base)
    }
    
    if(normalizeFunc == "quantile"){
        normalized_matrix = quantile_normalize(data_matrix)
    } else if(normalizeFunc == "medianCentering"){
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
