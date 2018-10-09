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

#' Median normalization of the data (per batch median)
#'
#' @name normalize
#'
#' @return
#' @export
#'
#' @examples
normalize_medians_batch <- function(data_long, sample_annotation = NULL,
                                    sample_id_col = 'FullRunName',
                                    batch_col = 'MS_batch.final',
                                    feature_id_col = 'peptide_group_label',
                                    measure_col = 'Intensity'){
  if (!(sample_id_col %in% names(data_long) & batch_col %in% names(data_long)) &
      !is.null(sample_annotation)){
    data_long = data_long %>% merge(sample_annotation)
  }
  df_normalized = data_long %>%
    group_by_at(vars(one_of(batch_col, feature_id_col))) %>%
    mutate(median_batch = median(UQ(sym(measure_col)), na.rm = T)) %>%
    ungroup() %>%
    mutate(median_global = median(UQ(sym(measure_col)), na.rm = T)) %>%
    mutate(diff = median_global - median_batch) %>%
    mutate(Intensity_normalized = UQ(sym(measure_col))+diff)

  return(df_normalized)
}

#' Median normalization of the data (global)
#'
#' @name normalize
#'
#' @return
#' @export
#'
#' @examples
normalize_medians_global <- function(data_long,
                                     sample_id_col = 'FullRunName',
                                    measure_col = 'Intensity'){
  df_normalized = data_long  %>%
    group_by_at(vars(one_of(sample_id_col))) %>%
    mutate(median_run = median(UQ(sym(measure_col)), na.rm = T)) %>%
    ungroup()
  df_normalized = df_normalized %>%
    mutate(median_global = median(UQ(sym(measure_col)), na.rm = T)) %>%
    mutate(diff = median_global - median_run) %>%
    mutate(Intensity_normalized = UQ(sym(measure_col))+diff)
  return(df_normalized)
}

#' normalize with the custom (continuous) fit
#'
#' @name normalize
#' @param sample_order_col column, determining the order of sample MS run, used
#'   as covariate to fit the non-linear fit
#' @param fit_func function to fit the (non)-linear trend
#' @param return_long whether the result should be the "long" data frame (as
#'   `df_long`) or "wide" (as `data_matrix`)
#' @param ... other parameters, usually those of the `fit_func`
#'
#' @return
#' @export
#'
#' @examples
#'
#' @seealso \code{\link{fit_nonlinear}}
normalize_custom_fit <- function(data_matrix, sample_annotation,
                                 batch_col = 'MS_batch.final',
                                 feature_id_col = 'peptide_group_label',
                                 sample_id_col = 'FullRunName',
                                 measure_col = 'Intensity',
                                 sample_order_col = 'order',
                                 fit_func = fit_nonlinear, ...){

  data_matrix = as.data.frame(data_matrix)
  data_matrix[[feature_id_col]] = rownames(data_matrix)
  #TODO: change to matrix_to_long
  df_long = data_matrix %>%
    melt(id.vars = feature_id_col)
  names(df_long) = c(feature_id_col, sample_id_col, measure_col)

  df_normalized = df_long %>%
    filter(!is.na(UQ(as.name(measure_col)))) %>% #filter(!is.na(Intensity))
    merge(sample_annotation) %>%
    arrange_(feature_id_col, sample_order_col) %>%
    group_by_at(vars(one_of(c(feature_id_col, batch_col)))) %>% #group_by(peptide_group_label, MS_batch.final) )
    filter(n() >2)%>%
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
    #TODO: try to get rid of rlang by the following expression:
    #mutate(Intensity_normalized = diff + UQ(sym(measure_col)))
    #if only the fitted data table is required (not recommended)
    fit_df = df_normalized %>% select(one_of(c('fit', feature_id_col,
                                               sample_id_col, batch_col)))

    casting_formula =  as.formula(paste(feature_id_col, sample_id_col,
                                        sep =  " ~ "))
    df_normalized = dcast(df_normalized, formula = casting_formula,
                          value.var = 'Intensity_normalized')
    df_normalized_matrix = as.matrix(df_normalized[,2:ncol(df_normalized)])
    rownames(df_normalized_matrix) = df_normalized[,1]

  return(list(data_matrix = df_normalized_matrix,
              fit_df = fit_df))
}


#' Standardized input-output ComBat normalization ComBat allows users to adjust
#' for batch effects in datasets where the batch covariate is known, using
#' methodology described in Johnson et al. 2007. It uses either parametric or
#' non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects.  Users are returned an expression matrix that has been corrected for
#' batch effects. The input data are assumed to be cleaned and normalized before
#' batch effect removal.
#'
#' @name normalize
#' @param par.prior
#'
#' @return `data_matrix`-size data matrix with batch-effect corrected by
#'   `ComBat`
#' @export
#'
#' @examples
correct_with_ComBat <- function(data_matrix, sample_annotation,
                                batch_col = 'MS_batch.final', par.prior = TRUE){
  batches = sample_annotation[[batch_col]]
  modCombat = model.matrix(~1, data = sample_annotation)
  corrected_proteome = sva::ComBat(dat = data_matrix, batch = batches,
                              mod = modCombat, par.prior = par.prior)
  return(corrected_proteome)
}
