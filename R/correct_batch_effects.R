#' @title Batch correction of normalized data
#' 
#' @description Batch correction of normalized data. Batch correction 
#' brings each feature in each batch to the comparable shape.
#' Currently the following batch correction functions are implemented:
#' \enumerate{
#'   \item Per-feature median centering: 
#'   \code{center_feature_batch_medians_df()}. 
#'   Median centering of the features (per batch median).
#'   \item correction with ComBat:  \code{correct_with_ComBat_df()}. 
#' Adjusts for discrete batch effects using ComBat. ComBat, described in 
#' Johnson et al. 2007. It uses either parametric or
#' non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects. Users are returned an expression matrix that has been corrected for
#' batch effects. The input data are assumed to be free of missing values
#' and normalized before batch effect removal. Please note that missing values 
#' are common in proteomics, which is why in some cases corrections like 
#' \code{center_peptide_batch_medians_df} are more appropriate.
#'   \item Continuous drift correction:  \code{adjust_batch_trend_df()}. 
#' Adjust batch signal trend with the custom (continuous) fit.
#' Should be followed by discrete corrections,
#' e.g. \code{center_feature_batch_medians_df()} or 
#' \code{correct_with_ComBat_df()}.
#' }
#' Alternatively, one can call the correction function with  
#' \code{correct_batch_effects_df()} wrapper. 
#' Batch correction method allows correction of 
#' continuous signal drift within batch (if required) and adjustment for
#' discrete difference across batches. 
#' 
#' 
#' @inheritParams proBatch
#' @param return_fit_df (logical) whether to return the \code{fit_df} from 
#' \code{adjust_batch_trend_dm} or only the data matrix
#' @param fit_func function to fit the (non)-linear trend
#' @param min_measurements the number of samples in a batch required for curve fitting.
#' @param par.prior use parametrical or non-parametrical prior 
#' @param continuous_func function to use for the fit (currently 
#' only \code{loess_regression} available); if order-associated fix is not 
#' required, should be \code{NULL}.
#' @param discrete_func function to use for adjustment of discrete batch effects
#' (\code{MedianCentering} or \code{ComBat}).
#' @param ... other parameters, usually of \code{adjust_batch_trend}, 
#' and \code{fit_func}.

#' 
#' @return the data in the same format as input (\code{data_matrix} or 
#' \code{df_long}).
#' For \code{df_long} the data frame stores the original values of 
#' \code{measure_col}
#' in another column called "preBatchCorr_[measure_col]", and the normalized 
#' values in \code{measure_col} column.
#' 
#' The function \code{adjust_batch_trend_dm()}, if \code{return_fit_df} is 
#' \code{TRUE} returns list of two items: 
#' \enumerate{
#'   \item \code{data_matrix}
#'   \item \code{fit_df}, used to examine the fitting curves
#'}

#' 
#' @examples 
#' 
#' #Median centering per feature per batch:
#' median_centered_df <- center_feature_batch_medians_df(
#' example_proteome, example_sample_annotation)
#' 
#' #Correct with ComBat: 
#' combat_corrected_df <- correct_with_ComBat_df(example_proteome, 
#' example_sample_annotation)
#' 
#' #Adjust the MS signal drift:
#' test_peptides = unique(example_proteome$peptide_group_label)[1:3]
#' test_peptide_filter = example_proteome$peptide_group_label %in% test_peptides
#' test_proteome = example_proteome[test_peptide_filter,]
#' adjusted_df <- adjust_batch_trend_df(test_proteome, 
#' example_sample_annotation, span = 0.7, 
#' min_measurements = 8)
#' plot_fit <- plot_with_fitting_curve(unique(adjusted_df$peptide_group_label), 
#' df_long = adjusted_df, measure_col = 'preTrendFit_Intensity',
#' fit_df = adjusted_df, sample_annotation = example_sample_annotation)
#' 
#' #Correct the data in one go:
#' batch_corrected_matrix <- correct_batch_effects_df(example_proteome, 
#' example_sample_annotation, 
#' continuous_func = 'loess_regression',
#' discrete_func = 'MedianCentering', 
#' batch_col = 'MS_batch',  
#' span = 0.7, min_measurements = 8)
#' 
#' @seealso \code{\link{fit_nonlinear}}
#' @name correct_batch_effects
NULL

#' 
#' @export
#' @rdname correct_batch_effects
#' 
center_feature_batch_medians_df <- function(df_long, sample_annotation = NULL,
                                         sample_id_col = 'FullRunName',
                                         batch_col = 'MS_batch',
                                         feature_id_col = 'peptide_group_label',
                                         measure_col = 'Intensity',
                                         keep_all = 'default',
                                         no_fit_imputed = TRUE,
                                         qual_col = NULL, 
                                         qual_value = NULL){
  
  original_cols = names(df_long)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, 
                                     order_col = NULL, 
                                     facet_col = NULL, merge = TRUE)
  
  old_measure_col = paste('preBatchCorr', measure_col, sep = '_')
  
  corrected_df = df_long %>%
    group_by_at(vars(one_of(batch_col, feature_id_col))) 
  
  if(is.null(qual_col) & no_fit_imputed){
    warning('imputed value flag column is NULL, changing no_fit_imputed to FALSE')
    no_fit_imputed = FALSE
  }
  
  if (no_fit_imputed){
    if(!is.null(qual_col) && !(qual_col %in% names(corrected_df))){
      stop("imputed value flag column (qual_col) is not in the data frame!")
    }
    temp_measure_col = paste('temp', measure_col, sep = '_')
    corrected_df = corrected_df %>%
      mutate(!!(sym(temp_measure_col)) := ifelse(!!sym(qual_col)== qual_value,
                                                 NA, !!(sym(measure_col)))) %>%
      mutate(median_batch = median(!!(sym(temp_measure_col)), na.rm = TRUE)) %>%
      ungroup() %>%
      group_by_at(vars(one_of(feature_id_col))) %>%
      mutate(median_global = median(!!(sym(temp_measure_col)), na.rm = TRUE))
  } else {
    if(!is.null(qual_col)){
      warning('imputed values are specified not to be used for median inference, 
however, 
              no_fit_imputed is FALSE')
    }
    corrected_df = corrected_df %>%
      mutate(median_batch = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
      ungroup() %>%
      group_by_at(vars(one_of(feature_id_col))) %>%
      mutate(median_global = median(!!(sym(measure_col)), na.rm = TRUE))
  }
  corrected_df = corrected_df %>%
    ungroup() %>%
    mutate(diff_medians = median_global - median_batch) %>%
    rename(!!(sym(old_measure_col)) := !!(sym(measure_col))) %>%
    mutate(!!(sym(measure_col)) := !!(sym(old_measure_col)) + diff_medians)
  
  if(!is.null(qual_col) && qual_col %in% names(corrected_df)){
    corrected_df = switch (keep_all,
                           all = corrected_df,
                           default = corrected_df %>%
                             dplyr::select(all_of(c(original_cols, old_measure_col, 
                                             qual_col, qual_value))),
                           minimal = corrected_df %>%
                             dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                             old_measure_col, qual_col, qual_value)))
    )
  } else {
    corrected_df = switch (keep_all,
                           all = corrected_df,
                           default = corrected_df %>%
                             dplyr::select(all_of(c(original_cols, old_measure_col))),
                           minimal = corrected_df %>%
                             dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                             old_measure_col)))
    )
  }
  
  
  return(corrected_df)
}

#' @export
#' @rdname correct_batch_effects
#' 
center_feature_batch_medians_dm <- function(data_matrix, sample_annotation,
                                            sample_id_col = 'FullRunName',
                                            batch_col = 'MS_batch',
                                            feature_id_col = 'peptide_group_label',
                                            measure_col = 'Intensity'){
  
  df_long = matrix_to_long(data_matrix, feature_id_col = feature_id_col,
                           measure_col = measure_col, 
                           sample_id_col = sample_id_col)
  
  corrected_df = center_feature_batch_medians_df(df_long, sample_annotation,
                                                 sample_id_col = sample_id_col,
                                                 batch_col = batch_col, 
                                                 feature_id_col = feature_id_col,
                                                 measure_col = measure_col)
  
  corrected_matrix = long_to_matrix(corrected_df, 
                                    sample_id_col = sample_id_col,
                                    measure_col = measure_col,
                                    feature_id_col = feature_id_col)
  return(corrected_matrix)
}

#' 
#' @export
#' @rdname correct_batch_effects
#' 
center_feature_batch_means_df <- function(df_long, sample_annotation = NULL,
                                            sample_id_col = 'FullRunName',
                                            batch_col = 'MS_batch',
                                            feature_id_col = 'peptide_group_label',
                                            measure_col = 'Intensity',
                                            keep_all = 'default',
                                            no_fit_imputed = TRUE,
                                            qual_col = NULL, 
                                            qual_value = NULL){
  
  original_cols = names(df_long)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, 
                                     order_col = NULL, 
                                     facet_col = NULL, merge = TRUE)
  
  old_measure_col = paste('preBatchCorr', measure_col, sep = '_')
  
  corrected_df = df_long %>%
    group_by_at(vars(one_of(batch_col, feature_id_col))) 
  
  if(is.null(qual_col) & no_fit_imputed){
    warning('imputed value flag column is NULL, changing no_fit_imputed to FALSE')
    no_fit_imputed = FALSE
  }
  
  if (no_fit_imputed){
    if(!is.null(qual_col) && !(qual_col %in% names(corrected_df))){
      stop("imputed value flag column (qual_col) is not in the data frame!")
    }
    temp_measure_col = paste('temp', measure_col, sep = '_')
    corrected_df = corrected_df %>%
      mutate(!!(sym(temp_measure_col)) := ifelse(!!sym(qual_col)== qual_value,
                                                 NA, !!(sym(measure_col)))) %>%
      mutate(mean_batch = mean(!!(sym(temp_measure_col)), na.rm = TRUE)) %>%
      ungroup() %>%
      group_by_at(vars(one_of(feature_id_col))) %>%
      mutate(mean_global = mean(!!(sym(temp_measure_col)), na.rm = TRUE))
  } else {
    if(!is.null(qual_col)){
      warning('imputed values are specified not to be used for mean inference, 
              however, 
              no_fit_imputed is FALSE')
    }
    corrected_df = corrected_df %>%
      mutate(mean_batch = mean(!!(sym(measure_col)), na.rm = TRUE)) %>%
      ungroup() %>%
      group_by_at(vars(one_of(feature_id_col))) %>%
      mutate(mean_global = mean(!!(sym(measure_col)), na.rm = TRUE))
    }
  corrected_df = corrected_df %>%
    ungroup() %>%
    mutate(diff_means = mean_global - mean_batch) %>%
    rename(!!(sym(old_measure_col)) := !!(sym(measure_col))) %>%
    mutate(!!(sym(measure_col)) := !!(sym(old_measure_col)) + diff_means)
  
  if(!is.null(qual_col) && qual_col %in% names(corrected_df)){
    corrected_df = switch (keep_all,
                           all = corrected_df,
                           default = corrected_df %>%
                             dplyr::select(all_of(c(original_cols, old_measure_col, 
                                             qual_col, qual_value))),
                           minimal = corrected_df %>%
                             dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                             old_measure_col, qual_col, qual_value)))
    )
  } else {
    corrected_df = switch (keep_all,
                           all = corrected_df,
                           default = corrected_df %>%
                             dplyr::select(all_of(c(original_cols, old_measure_col))),
                           minimal = corrected_df %>%
                             dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                             old_measure_col)))
    )
  }
  
  
  return(corrected_df)
}

#' @export
#' @rdname correct_batch_effects
#' 
center_feature_batch_means_dm <- function(data_matrix, sample_annotation,
                                            sample_id_col = 'FullRunName',
                                            batch_col = 'MS_batch',
                                            feature_id_col = 'peptide_group_label',
                                            measure_col = 'Intensity'){
  
  df_long = matrix_to_long(data_matrix, feature_id_col = feature_id_col,
                           measure_col = measure_col, 
                           sample_id_col = sample_id_col)
  
  corrected_df = center_feature_batch_means_df(df_long, sample_annotation,
                                                 sample_id_col = sample_id_col,
                                                 batch_col = batch_col, 
                                                 feature_id_col = feature_id_col,
                                                 measure_col = measure_col)
  corrected_matrix = long_to_matrix(corrected_df, 
                                    sample_id_col = sample_id_col,
                                    measure_col = measure_col,
                                    feature_id_col = feature_id_col)
  
  return(corrected_matrix)
}

#' 
#' @export
#' @rdname correct_batch_effects
#'
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{plot_with_fitting_curve}}
adjust_batch_trend_df <- function(df_long, sample_annotation = NULL,
                               batch_col = 'MS_batch',
                               feature_id_col = 'peptide_group_label',
                               sample_id_col = 'FullRunName',
                               measure_col = 'Intensity',
                               order_col = 'order',
                               keep_all = 'default',
                               fit_func = 'loess_regression', 
                               no_fit_imputed = TRUE,
                               qual_col = NULL, 
                               qual_value = NULL,
                               min_measurements = 8, ...){
  
  original_cols = names(df_long)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, 
                                     facet_col = NULL, merge = TRUE)
  
  if(!is.null(qual_col) & no_fit_imputed){
    if(!(qual_col %in% names(df_long))){
      stop("imputed value flag column is not in the data frame!")
    }
  
  } else{
    if(is.null(qual_col) & no_fit_imputed){
      warning('imputed value flag column is NULL, changing no_fit_imputed to FALSE')
      no_fit_imputed = FALSE
    }
  }
  
  if (!is.null(batch_col)){
    if (!(batch_col %in% union(names(df_long),names(sample_annotation)))){
      stop('Batch column is neither in sample annotation nor in data matrix')
    }
    sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
    corrected_df = df_long %>%
      #filter(!is.na(!!(sym(measure_col)))) %>% #filter(!is.na(Intensity))
      group_nest(!!!syms(c(feature_id_col, batch_col))) %>%  
      mutate(fit = pmap(list(df_feature_batch = data, 
                             feature_id = !!sym(feature_id_col), 
                             batch_id = as.character(!!sym(batch_col))),
                        fit_nonlinear, 
                        measure_col = measure_col,order_col = order_col, 
                        fit_func = fit_func, 
                        no_fit_imputed = no_fit_imputed,
                        qual_col = qual_col, qual_value = qual_value,
                        min_measurements = min_measurements, ...))
  } else {
    corrected_df = df_long %>%
      #filter(!is.na(!!(sym(measure_col)))) %>% #filter(!is.na(Intensity))
      nest(!!sym(feature_id_col)) %>%
      mutate(fit = map(data, fit_nonlinear, 
                        measure_col = measure_col,order_col = order_col, 
                        fit_func = fit_func,
                       min_measurements = min_measurements, ...))
  }
  
  old_measure_col = paste('preTrendFit', measure_col, sep = '_')
  
  corrected_df = corrected_df %>%
    unnest() %>%
    group_by(!!!syms(c(feature_id_col, batch_col))) %>%
    mutate(mean_fit = mean(fit, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(diff_fit = mean_fit - fit) %>%
    mutate(diff.na = ifelse(is.na(diff_fit), 0, diff_fit)) %>%
    rename(!!(old_measure_col) := !!(sym(measure_col))) %>%
    #TODO: make the conditional shift: if we keep the requants, 
    #then we can fit the curve without, but shift them nevertheless
    mutate(!!(sym(measure_col)) := !!sym('diff.na') + !!sym(old_measure_col))
  
  if(!is.null(qual_col) && qual_col %in% names(corrected_df)){
    corrected_df = switch (keep_all,
                           all = corrected_df,
                           default = corrected_df %>%
                             dplyr::select(all_of(c(original_cols, old_measure_col, 
                                             'fit', batch_col, qual_col, qual_value))),
                           minimal = corrected_df %>%
                             dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                             old_measure_col,'fit', batch_col, qual_col, qual_value)))
    )
  } else {
  
  corrected_df = switch (keep_all,
    all = corrected_df ,
    default = corrected_df %>%
      dplyr::select(all_of(c(original_cols, old_measure_col, 'fit'))),
    minimal = corrected_df %>%
      dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                      old_measure_col,'fit')))
  )
  }
  
  return(corrected_df)
  }

#' 
#' @export
#' @rdname correct_batch_effects
#'
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{plot_with_fitting_curve}}
adjust_batch_trend_dm <- function(data_matrix, sample_annotation,
                                  batch_col = 'MS_batch',
                                  feature_id_col = 'peptide_group_label',
                                  sample_id_col = 'FullRunName',
                                  measure_col = 'Intensity',
                                  order_col = 'order',
                                  fit_func = 'loess_regression', 
                                  return_fit_df = TRUE,
                                  min_measurements = 8, ...){
  df_long = matrix_to_long(data_matrix, feature_id_col = feature_id_col,
                           measure_col = measure_col, 
                           sample_id_col = sample_id_col)
  
  corrected_data = adjust_batch_trend_df(df_long, sample_annotation,
                                       sample_id_col = sample_id_col,
                                       batch_col = batch_col, 
                                       feature_id_col = feature_id_col,
                                       measure_col = measure_col)
  
  corrected_df = corrected_data$corrected_df
  corrected_dm = long_to_matrix(corrected_df, feature_id_col = feature_id_col,
                                measure_col = measure_col, 
                                sample_id_col = sample_id_col)
  if (return_fit_df){
    
    fit_df = corrected_data$fit_df
    return(list(corrected_dm = corrected_dm,
                fit_df = fit_df))
  } else {
    return(corrected_dm)
  }
  
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
                                   par.prior = TRUE,
                                   no_fit_imputed = TRUE,
                                   qual_col = NULL, 
                                   qual_value = NULL,
                                   keep_all = 'default'){
  
  original_cols = names(df_long)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, 
                                     facet_col = NULL, merge = FALSE)
  
  #TODO: no_impute_values fix
  if(!is.null(qual_col) & no_fit_imputed){
    if(!(qual_col %in% names(df_long))){
      stop("imputed value flag column is not in the data frame!")
    }
    data_matrix = long_to_matrix(df_long, feature_id_col = feature_id_col,
                                 measure_col = measure_col, 
                                 sample_id_col = sample_id_col, 
                                 qual_col = qual_col, qual_value = qual_value)
  } else {
    data_matrix = long_to_matrix(df_long, feature_id_col = feature_id_col,
                                 measure_col = measure_col, 
                                 sample_id_col = sample_id_col, 
                                 qual_col = NULL)
  }
  
  
  corrected_matrix = run_ComBat_core(sample_annotation, batch_col, data_matrix, 
                                     par.prior)
  
  corrected_df = matrix_to_long(corrected_matrix, 
                                feature_id_col = feature_id_col,
                                measure_col = measure_col, 
                                sample_id_col = sample_id_col)
  
  old_measure_col = paste('preBatchCorr', measure_col, sep = '_')
  
  df_long = df_long %>% 
    rename(!!(old_measure_col) := !!(sym(measure_col))) 
  
  corrected_df = corrected_df %>%
    merge(df_long, by = c(feature_id_col, sample_id_col))
  
  corrected_df = switch (keep_all,
                         all = corrected_df ,
                         default = corrected_df %>%
                           dplyr::select(all_of(c(original_cols, old_measure_col))),
                         minimal = corrected_df %>%
                           dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                           old_measure_col))))
  
  return(corrected_df)
}

#' 
#' @export
#' @rdname correct_batch_effects
#'
correct_with_ComBat_dm <- function(data_matrix, sample_annotation = NULL,
                                   feature_id_col = 'peptide_group_label',
                                   measure_col = 'Intensity',
                                   sample_id_col = 'FullRunName',
                                   batch_col = 'MS_batch', 
                                   par.prior = TRUE){
  
  df_long = matrix_to_long(data_matrix, feature_id_col = feature_id_col,
                           measure_col = measure_col, 
                           sample_id_col = sample_id_col)
  
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, 
                                     facet_col = NULL, merge = FALSE)
  
  data_matrix = long_to_matrix(df_long, feature_id_col = feature_id_col,
                               measure_col = measure_col, 
                               sample_id_col = sample_id_col)
  
  corrected_matrix = run_ComBat_core(sample_annotation, batch_col, data_matrix, 
                                     par.prior)
  
  return(corrected_matrix)
}

run_ComBat_core <- function(sample_annotation, batch_col, data_matrix, 
                            par.prior) {
  #TODO: program for the case of multiple batch factors - "SuperBatch"
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
correct_batch_effects_df <- function(df_long, sample_annotation, 
                                  continuous_func = NULL, 
                                  discrete_func = c("MedianCentering","MeanCentering","ComBat"),
                                  batch_col = 'MS_batch',  
                                  feature_id_col = 'peptide_group_label', 
                                  sample_id_col = 'FullRunName',
                                  measure_col = 'Intensity',  
                                  order_col = 'order', 
                                  keep_all = 'default',
                                  no_fit_imputed = TRUE,
                                  qual_col = NULL, 
                                  qual_value = NULL,
                                  min_measurements = 8, ...){
  
  discrete_func <- match.arg(discrete_func)
  
  sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
  
  original_cols = names(df_long)
  
  if(!is.null(continuous_func)){
    df_long = adjust_batch_trend_df(df_long = df_long, 
                                  sample_annotation = sample_annotation,
                                  batch_col = batch_col,
                                  feature_id_col = feature_id_col,
                                  sample_id_col = sample_id_col,
                                  measure_col = measure_col,
                                  order_col = order_col, 
                                  keep_all = keep_all,
                                  no_fit_imputed = no_fit_imputed,
                                  qual_col = qual_col, 
                                  qual_value = qual_value,
                                  fit_func = continuous_func, 
                                  min_measurements = min_measurements, ...)
    
  }
  
  #TODO: re-implement in a functional programming
  if(discrete_func == 'MedianCentering'){
    corrected_df = center_feature_batch_medians_df(df_long = df_long, 
                                               sample_annotation = sample_annotation,
                                               sample_id_col = sample_id_col,
                                               batch_col = batch_col,
                                               feature_id_col = feature_id_col,
                                               measure_col = measure_col,
                                               no_fit_imputed = no_fit_imputed,
                                               qual_col = qual_col, 
                                               qual_value = qual_value)
  } 
  
  if(discrete_func == 'MeanCentering'){
    corrected_df = center_feature_batch_means_df(df_long = df_long, 
                                                   sample_annotation = sample_annotation,
                                                   sample_id_col = sample_id_col,
                                                   batch_col = batch_col,
                                                   feature_id_col = feature_id_col,
                                                   measure_col = measure_col,
                                                   no_fit_imputed = no_fit_imputed,
                                                   qual_col = qual_col, 
                                                   qual_value = qual_value)
  } 
  
  if(discrete_func == 'ComBat'){
    #TODO: pick up the functional programming argument borrowing
    corrected_df = correct_with_ComBat_df(df_long = df_long, 
                                          sample_annotation = sample_annotation,
                                          feature_id_col = feature_id_col, 
                                          measure_col = measure_col,
                                          sample_id_col = sample_id_col,
                                          batch_col = batch_col, 
                                          par.prior = TRUE)
    #TODO: fix for "no_fit_imputed" cases
  }
  
  old_measure_col = paste('preBatchCorr', measure_col, sep = '_')
  if(!is.null(continuous_func)){
    preFit_measure_col = paste('preTrendFit', measure_col, sep = '_')
    corrected_df = switch (keep_all,
                           all = corrected_df ,
                           default = corrected_df %>%
                             dplyr::select(all_of(c(original_cols, old_measure_col, preFit_measure_col, 'fit'))),
                           minimal = corrected_df %>%
                             dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                             old_measure_col)))
    )} else {
      corrected_df = switch (keep_all,
                             all = corrected_df ,
                             default = corrected_df %>%
                               dplyr::select(all_of(c(original_cols, old_measure_col))),
                             minimal = corrected_df %>%
                               dplyr::select(all_of(c(sample_id_col, feature_id_col, measure_col, 
                                               old_measure_col))))
  }
  
  return(corrected_df)
}

#' 
#' @export
#' @rdname correct_batch_effects
#' 
correct_batch_effects_dm <- function(data_matrix, sample_annotation, 
                                     continuous_func = NULL, 
                                     discrete_func = c("MedianCentering", 
                                                       "ComBat"),
                                     batch_col = 'MS_batch',  
                                     feature_id_col = 'peptide_group_label', 
                                     sample_id_col = 'FullRunName',
                                     measure_col = 'Intensity',  
                                     order_col = 'order', 
                                     min_measurements = 8, ...){
  
  df_long = matrix_to_long(data_matrix, feature_id_col = feature_id_col,
                           measure_col = measure_col, 
                           sample_id_col = sample_id_col)
  corrected_df = correct_batch_effects_df(df_long, sample_annotation, 
                                          continuous_func = continuous_func, 
                                          discrete_func = discrete_func,
                                          batch_col = batch_col,  
                                          feature_id_col = feature_id_col, 
                                          sample_id_col = sample_id_col,
                                          measure_col = measure_col, 
                                          order_col = order_col,
                                          min_measurements = min_measurements, 
                                          no_fit_imputed = TRUE,
                                          qual_col = NULL, 
                                          qual_value = NULL, 
                                          keep_all = FALSE,...)
  corrected_matrix = long_to_matrix(corrected_df, 
                                    sample_id_col = sample_id_col,
                                    measure_col = measure_col,
                                    feature_id_col = feature_id_col)
  
  return(corrected_matrix)
}
