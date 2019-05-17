#' Batch correction of normalized data
#' 
#' @description Batch correction of normalized data. Batch correction 
#' brings each feature in each batch to the comparable shape.
#' Currently the following batch correction functions are implemented:
#' \enumerate{
#'   \item Per-feature median centering: `center_feature_batch_medians()`. Median centering of the features (per batch median)
#'   \item correction with ComBat: `correct_with_ComBat()`. Adjusts for discrete batch effects using ComBat. 
#'   Standardized input-output ComBat normalization ComBat allows users to adjust
#' for batch effects in datasets where the batch covariate is known, using
#' methodology described in Johnson et al. 2007. It uses either parametric or
#' non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects. Users are returned an expression matrix that has been corrected for
#' batch effects. The input data are assumed to be free of missing values
#'  and normalized before batch effect removal. Please note that missing values 
#'  are common in proteomics, which is why in some cases corrections like 
#'  \code{center_peptide_batch_medians} are more appropriate.
#'   \item Continuous drift correction: `adjust_batch_trend()`. Adjust batch signal 
#'   trend with the custom (continuous) fit. Should be followed by discrete corrections,
#'   e.g. `center_feature_batch_medians()` or `correct_with_ComBat()`.
#' }
#' Alternatively, one can call the correction function with `correct_batch_effects()` wrapper. #' Batch correction method allows correction of 
#' continuous signal drift within batch (if required) and adjustment for
#' discrete difference across batches. 
#' 
#' 
#' @inheritParams proBatch
#' @param fit_func function to fit the (non)-linear trend
#' @param abs_threshold the absolute threshold (number of samples in a batch) to filter data for curve fitting 
#' @param pct_threshold the percentage threshold (fraction of samples in a batch)  to filter data for curve fitting 
#' @param par.prior whether parametrical or non-parametrical prior should be used
#' @param continuous_func function to use for the fit (currently 
#' only \code{loess_regression} available); if order-associated fix is not required, should be `NULL`
#' @param discrete_func function to use for adjustment of discrete batch effects (\code{MedianCentering} or \code{ComBat})
#' @param ... other parameters, usually of \code{adjust_batch_trend}, and \code{fit_func}

#' 
#' @return the data in the same format as input (\code{data_matrix} or \code{df_long}).
#' For \code{df_long} the data frame stores the original values of \code{measure_col}
#' in another column called "old_intensity" if "intensity", and the normalized values
#' in \code{measure_col} column.
#' 
#' The function `adjust_batch_trend()` returns list of two items: 
#' \enumerate{
#'   \item \code{data_matrix}
#'   \item \code{fit_df}, used to examine the fitting curves
#'}

#' 
#' @examples 
#' 
#' #Median centering per feature per batch:
#' median_centered_df <- center_feature_batch_medians(
#' example_proteome, example_sample_annotation)
#' 
#' #Correct with ComBat: 
#' combat_corrected_df <- correct_with_ComBat(example_proteome, example_sample_annotation)
#' 
#' #Adjust the MS signal drift:
#' adjusted_data <- adjust_batch_trend(example_proteome[example_proteome$peptide_group_label %in% unique(example_proteome$peptide_group_label)[1:3],], 
#' example_sample_annotation, span = 0.7, 
#' abs_threshold = 5, pct_threshold = 0.20)
#' fit_df <- adjusted_data$fit_df
#' adjusted_df <- adjusted_data$corrected_df
#' plot_fit <- plot_with_fitting_curve(unique(adjusted_df$peptide_group_label), 
#' df_long = adjusted_df, fit_df = fit_df, 
#' sample_annotation = example_sample_annotation)
#' 
#' #Correct the data in one go:
#' batch_corrected_matrix <- correct_batch_effects(example_proteome, example_sample_annotation, 
#' continuos_func = 'loess_regression',
#' discrete_func = 'MedianCentering', 
#' batch_col = 'MS_batch',  
#' span = 0.7,
#' abs_threshold = 5, pct_threshold = 0.20)
#' 
#' @seealso \code{\link{fit_nonlinear}}
#' @name correct_batch_effects
NULL

#' 
#' @export
#' @rdname correct_batch_effects
#' 
center_feature_batch_medians <- function(df_long, sample_annotation = NULL,
                                         sample_id_col = 'FullRunName',
                                         batch_col = 'MS_batch',
                                         feature_id_col = 'peptide_group_label',
                                         measure_col = 'Intensity'){
  
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, facet_col = NULL)
  
  old_measure_col <- paste('old', measure_col, sep = '_')
  corrected_df = df_long %>%
    group_by_at(vars(one_of(batch_col, feature_id_col))) %>%
    mutate(median_batch = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup() %>%
    group_by_at(vars(one_of(feature_id_col))) %>%
    mutate(median_global = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(diff = median_global - median_batch) %>%
    rename(!!(old_measure_col) := !!(sym(measure_col))) %>%
    mutate(!!(sym(measure_col)) := !!(sym(old_measure_col))+diff)
  
  return(corrected_df)
}

#' 
#' @export
#' @rdname correct_batch_effects
#'
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{plot_with_fitting_curve}}
adjust_batch_trend <- function(df_long, sample_annotation = NULL,
                                 batch_col = 'MS_batch',
                                 feature_id_col = 'peptide_group_label',
                                 sample_id_col = 'FullRunName',
                                 measure_col = 'Intensity',
                                 order_col = 'order',
                                 fit_func = 'loess_regression', 
                                 abs_threshold = 5, pct_threshold = 0.20, ...){
  
  sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
  
  sample_annotation = sample_annotation %>%
    group_by(!!sym(batch_col)) %>%
    mutate(batch_total = n()) %>%
    ungroup() %>%
    mutate(!!batch_col := as.character(!!sym(batch_col)))
  
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, batch_col, order_col, facet_col = NULL)
  
  corrected_df = df_long %>%
    filter(!is.na(!!(sym(measure_col)))) %>% #filter(!is.na(Intensity))
    group_nest(!!!syms(c(feature_id_col, batch_col, "batch_total"))) %>%  
    mutate(fit = pmap(list(df_feature_batch = data,  batch_size = batch_total, 
                                  feature_id = !!sym(feature_id_col)), batch_id = !!sym(batch_col),
                             fit_nonlinear, 
                             measure_col = measure_col,order_col = order_col, 
                             fit_func = fit_func, 
                             abs_threshold = abs_threshold,pct_threshold = pct_threshold, ...)) %>%
    unnest() %>%
    group_by(!!!syms(c(feature_id_col, batch_col))) %>%
    mutate(mean_fit = mean(fit, na.rm = T)) %>%
    ungroup() %>%
    mutate(diff = mean_fit - fit) %>%
    mutate(diff.na = ifelse(is.na(diff), 0, diff)) %>%
    mutate_(Intensity_normalized = interp(~`+`(x, y),
                                          x = as.name('diff.na'),
                                          y = as.name(measure_col))) %>%
    rename(!!(paste('old', measure_col, sep = '_')) := !!(sym(measure_col))) %>%
    rename(!!(sym(measure_col)) := Intensity_normalized)
  
  fit_df = corrected_df %>% dplyr::select(one_of(c('fit', 'diff', 'diff.na', feature_id_col,
                                                    sample_id_col, measure_col, batch_col)))
  
  corrected_df = corrected_df %>%
    select(c(sample_id_col, feature_id_col, measure_col, paste('old', measure_col, sep = '_')))
  
  return(list(corrected_df = corrected_df,
              fit_df = fit_df))
  }

#' 
#' @export
#' @rdname correct_batch_effects
#'
correct_with_ComBat_df <- function(df_long, sample_annotation = NULL, 
                                feature_id_col = 'peptide_group_label',
                                measure_col = 'Intensity',
                                sample_id_col = 'FullRunName',
                                batch_col = 'MS_batch', 
                                par.prior = TRUE){
  
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, facet_col = NULL)
  
  data_matrix = long_to_matrix(df_long, feature_id_col = feature_id_col,
                               measure_col = measure_col, sample_id_col = sample_id_col)
  
  corrected_matrix = run_ComBat_core(sample_annotation, batch_col, data_matrix, par.prior)
  
  corrected_df = matrix_to_long(corrected_matrix, feature_id_col = feature_id_col,
                           measure_col = measure_col, sample_id_col = sample_id_col)
  
  return(corrected_df)
}

correct_with_ComBat_matrix <- function(data_matrix, sample_annotation = NULL, 
                                feature_id_col = 'peptide_group_label',
                                measure_col = 'Intensity',
                                sample_id_col = 'FullRunName',
                                batch_col = 'MS_batch', 
                                par.prior = TRUE){
  
  df_long = matrix_to_long(data_matrix, feature_id_col = feature_id_col,
                               measure_col = measure_col, sample_id_col = sample_id_col)
  
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, facet_col = NULL)
  
  data_matrix = long_to_matrix(df_long, feature_id_col = feature_id_col,
                               measure_col = measure_col, sample_id_col = sample_id_col)
  
  corrected_matrix = run_ComBat_core(sample_annotation, batch_col, data_matrix, par.prior)
  
  return(corrected_matrix)
}

run_ComBat_core <- function(sample_annotation, batch_col, data_matrix, par.prior) {
  batches = sample_annotation[[batch_col]]
  modCombat = model.matrix(~1, data = sample_annotation)
  corrected_matrix = ComBat(dat = data_matrix, batch = batches,
                            mod = modCombat, par.prior = par.prior)
  return(corrected_matrix)
}

#' 
#' @export
#' @rdname correct_batch_effects
#' 
correct_batch_effects <- function(df_long, sample_annotation, 
                                  continuous_func = NULL, 
                                  discrete_func = c("MedianCentering", "ComBat"), 
                                  batch_col = 'MS_batch',  
                                  feature_id_col = 'peptide_group_label', 
                                  sample_id_col = 'FullRunName',
                                  measure_col = 'Intensity',  
                                  sample_order_col = 'order', 
                                  abs_threshold = 5, pct_threshold = 0.20, ...){
  
  discrete_func <- match.arg(discrete_func)
  
  sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
  
  if(!is.null(continuous_func)){
    fit_list = adjust_batch_trend(df_long, 
                                  sample_annotation = sample_annotation,
                                  batch_col = batch_col,
                                  feature_id_col = feature_id_col,
                                  sample_id_col = sample_id_col,
                                  measure_col = measure_col,
                                  sample_order_col = sample_order_col,
                                  fit_func = continuous_func, 
                                  abs_threshold = abs_threshold, 
                                  pct_threshold = pct_threshold, ...)
    fit_matrix = fit_list$data_matrix
    fit_long = matrix_to_long(fit_matrix, feature_id_col = feature_id_col,
                              measure_col = measure_col, sample_id_col = sample_id_col)
    df_long = fit_long
  }
  
  #TODO: re-implement in a functional programming
  if(discrete_func == 'MedianCentering'){
    corrected_df = center_feature_batch_medians(df_long = df_long, 
                                               sample_annotation = sample_annotation,
                                               sample_id_col = sample_id_col,
                                               batch_col = batch_col,
                                               feature_id_col = feature_id_col,
                                               measure_col = measure_col)
  }
  
  if(discrete_func == 'ComBat'){
    #TODO: pick up the functional programming argument borrowing
    corrected_df = correct_with_ComBat_df(df_long = df_long, 
                                            sample_annotation = sample_annotation,
                                            feature_id_col = feature_id_col, 
                                            measure_col = measure_col,
                                            sample_id_col = sample_id_col,
                                            batch_col = batch_col, par.prior = TRUE)
  }
  
  return(corrected_df)
}

