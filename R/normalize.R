


#' Quantile normalization of the data, ensuring that the row and column names are retained
#'
#' @param data_matrix log transformed data matrix (features in rows and samples in columns)
#'
#' @return
#' @export
#' @import preprocessCore
#'
#' @examples
quantile_normalize <- function(data_matrix){
  q_norm_proteome = normalize.quantiles(data_matrix)
  colnames(q_norm_proteome) = colnames(data_matrix)
  rownames(q_norm_proteome) = rownames(data_matrix)
  return(q_norm_proteome)
}

#' Median normalization of the data
#'
#' @param sample_annotation
#' @param batch_column
#' @param measure_column
#' @param data_matrix
#'
#' @return
#' @export
#' @import tidyverse
#' @import lazyeval
#'
#' @examples
median_normalization <- function(data_matrix, sample_annotation,
                                 batch_column, measure_column)
  {
  data_matrix = data_matrix %>% merge(sample_annotation) %>%
    group_by_at(vars(one_of(batch_column))) %>%
    mutate_(median = interp(~median(measurement),
            measurement = as.name(measure_column))) %>%
    ungroup() %>%
    mutate_(median_overall = interp(~median(measurement),
            measurement = as.name(measure_column))) %>%
    mutate(diff = median - median_overall) %>%
    mutate()
}

#' normalize with the custom (continuous) fit
#'
#' @param data_matrix
#' @param sample_annotation
#' @param batch_col
#' @param feature_id_col
#' @param sample_id_column
#' @param measure_col
#' @param sample_order_col
#' @param fit_func
#' @param return_long
#' @param ... other parameters, usually those of the `fit_func`
#'
#' @return
#' @export
#' @import dplyr
#' @importFrom tidyr nest
#' @importFrom tidyr unnest
#' @import reshape2
#' @import lazyeval
#' @importFrom  purrr map
#'
#' @examples
normalize_custom_fit <- function(data_matrix, sample_annotation, batch_col,
                                 feature_id_col, sample_id_column, measure_col,
                                 sample_order_col, fit_func, return_long = F, ...){
  data_matrix = as.data.frame(data_matrix)
  data_matrix[[feature_id_col]] = rownames(data_matrix)
  df_long = data_matrix %>%
    melt(id.vars = feature_id_col)
  names(df_long) = c(feature_id_col, sample_id_column, measure_col)

  df_normalized = df_long %>%
    merge(sample_annotation) %>%
    arrange_(feature_id_col, sample_order_col) %>%
    group_by_at(vars(one_of(c(feature_id_col, batch_col)))) %>%
    nest() %>%
    mutate(fit = map(data, fit_func, response.var = measure_col,
                     expl.var = sample_order_col, ...)) %>%
    unnest() %>%
    #change the fit to the corrected data
    group_by_at(vars(one_of(c(feature_id_col, batch_col)))) %>%
    mutate(mean_fit = mean(fit)) %>%
    mutate(diff = mean_fit - fit) %>%
    mutate_(Intensity_normalized = interp(~`+`(x, y),
                                          x = as.name('diff'),
                                          y = as.name(measure_col)))
    #if only the fitted data table is required (not recommended)
    if(!return_long){
      casting_formula =  as.formula(paste(feature_id_col, sample_id_column,
                                          sep =  " ~ "))
      df_normalized = dcast(df_normalized, formula = casting_formula,
                            value.var = 'Intensity_normalized')
                     }
  return(df_normalized)
}


#' Standardized input-output ComBat normalization
#'
#' @param sample_annotation
#' @param batch_column
#' @param par.prior
#' @param data_matrix
#'
#' @return
#' @export
#' @import sva
#'
#' @examples
correct_with_ComBat <- function(data_matrix, sample_annotation,
                                batch_column = 'MS_batch.final', par.prior = TRUE){
  batches = sample_annotation[[batch_column]]
  modCombat = model.matrix(~1, data = sample_annotation)
  corrected_proteome = ComBat(dat = data_matrix, batch = batches,
                              mod = modCombat, par.prior = par.prior)
  return(corrected_proteome)
}
