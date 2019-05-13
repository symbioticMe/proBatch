#' Median centering of the features (per batch median)
#' 
#' @inheritParams proBatch
#' 
#' @return \code{df_long}-size long format data with batch-effect corrected with
#'   per-feature batch median centering in \code{Intensity_normalized} column
#'   
#' @examples 
#' median_centered_proteome <- center_peptide_batchmedians(
#' example_proteome, example_sample_annotation)
#' 
#' @export
#'
#'@seealso \code{\link{correct_batch_effects}}
center_peptide_batch_medians <- function(df_long, sample_annotation = NULL,
                                         sample_id_col = 'FullRunName',
                                         batch_col = 'MS_batch',
                                         feature_id_col = 'peptide_group_label',
                                         measure_col = 'Intensity'){
  
  if(!setequal(unique(sample_annotation[[sample_id_col]]), 
               unique(df_long[[sample_id_col]]))){
    warning('Sample IDs in sample annotation not 
                consistent with samples in input data.')}
  
  if (!(sample_id_col %in% names(df_long) & batch_col %in% names(df_long)) &
      !is.null(sample_annotation)){
    df_long = df_long %>% merge(sample_annotation, by = sample_id_col)
  }
  df_normalized = df_long %>%
    group_by_at(vars(one_of(batch_col, feature_id_col))) %>%
    mutate(median_batch = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup() %>%
    group_by_at(vars(one_of(feature_id_col))) %>%
    mutate(median_global = median(!!(sym(measure_col)), na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(diff = median_global - median_batch) %>%
    mutate(Intensity_normalized = !!(sym(measure_col))+diff)
  
  return(df_normalized)
}


#' adjust batch signal trend with the custom (continuous) fit
#'
#' @inheritParams proBatch
#' @param fit_func function to fit the (non)-linear trend
#' @param abs_threshold the absolute threshold (number of samples in a batch) to filter data for curve fitting 
#' @param pct_threshold the percentage threshold (fraction of samples in a batch)  to filter data for curve fitting 
#' @param ... other parameters, usually those of the \code{fit_func}
#'
#' @return list of two items: 1) `data_matrix`, adjusted with continuous fit; 
#' 2) fit_df, used to examine the fitting curves
#' @examples 
#' adjusted_data <- adjust_batch_trend(example_proteome %>% 
#' filter(peptide_group_label %in% unique(example_proteome$peptide_group_label)[1:3]), 
#' example_sample_annotation, span = 0.7, 
#' abs_threshold = 5, pct_threshold = 0.20)
#' fit_df <- adjusted_data$fit_df
#' adjusted_data_matrix <- adjusted_data$data_matrix
#' adjusted_df <- matrix_to_long(adjusted_data_matrix)
#' plot_fit <- plot_with_fitting_curve(unique(adjusted_df$peptide_group_label), 
#' df_long = adjusted_df, fit_df = fit_df, 
#' sample_annotation = example_sample_annotation)
#' 
#' @export
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
  
  #TODO: substitute with "check for sample consistency"
  sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
  sampleNames <- unique(df_long[[sample_id_col]])
  s_a <- sample_annotation[[sample_id_col]]
  all <- union(sampleNames, s_a)
  non_matched <- all[!all %in% intersect(sampleNames, s_a)]
  if(length(non_matched)!=0){
    warning("Sample ID in data matrix and sample annotation don't match. 
            Non-matching elements are removed for analysis")
  }
  
  sample_annotation = sample_annotation %>%
    filter(!!(as.name(sample_id_col)) %in% sampleNames) %>%
    arrange(match(!!(as.name(sample_id_col)), sampleNames)) %>%
    droplevels()
  
  sample_annotation = sample_annotation %>%
    group_by(!!sym(batch_col)) %>%
    mutate(batch_total = n()) %>%
    ungroup() %>%
    mutate(!!batch_col := as.character(!!sym(batch_col)))
  
  df_normalized = df_long %>%
    filter(!is.na(!!(sym(measure_col)))) %>% #filter(!is.na(Intensity))
    merge(sample_annotation, by = sample_id_col) %>%
    group_nest(!!!syms(c(feature_id_col, batch_col, "batch_total"))) %>%  
    mutate(fit = pmap(list(df_feature_batch = data,  batch_size = batch_total, 
                                  feature_id = !!sym(feature_id_col)), batch_id = !!sym(batch_col),
                             fit_nonlinear, 
                             measure_col = measure_col,order_col = order_col, 
                             fit_func = fit_func, 
                             abs_threshold = abs_threshold,pct_threshold = pct_threshold, ...)) %>%
    unnest() %>%
    group_by(!!!syms(c(feature_id_col, batch_col))) %>%
    mutate(mean_fit = mean(fit), na.rm = T) %>%
    ungroup() %>%
    mutate(diff = mean_fit - fit) %>%
    mutate(diff.na = ifelse(is.na(diff), 0, diff)) %>%
    mutate_(Intensity_normalized = interp(~`+`(x, y),
                                          x = as.name('diff.na'),
                                          y = as.name(measure_col)))
  
  fit_df = df_normalized %>% dplyr::select(one_of(c('fit', 'diff', 'diff.na', feature_id_col,
                                                    sample_id_col, measure_col)))
  
  casting_formula =  as.formula(paste(feature_id_col, sample_id_col,
                                      sep =  " ~ "))
  df_normalized = dcast(df_normalized, formula = casting_formula,
                        value.var = 'Intensity_normalized')
  df_normalized_matrix = as.matrix(df_normalized[,2:ncol(df_normalized)])
  rownames(df_normalized_matrix) = df_normalized[,1]
  
  return(list(data_matrix = df_normalized_matrix,
              fit_df = fit_df))
  }


#' Adjusts for discrete batch effects using ComBat
#'
#' @description Standardized input-output ComBat normalization ComBat allows users to adjust
#' for batch effects in datasets where the batch covariate is known, using
#' methodology described in Johnson et al. 2007. It uses either parametric or
#' non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects. Users are returned an expression matrix that has been corrected for
#' batch effects. The input data are assumed to be free of missing values
#'  and normalized before batch effect removal. Please note that missing values 
#'  are common in proteomics, which is why in some cases corrections like 
#'  \code{center_peptide_batch_medians} are more appropriate.
#'
#' @inheritParams proBatch
#' @param par.prior whether parametrical or non-parametrical prior should be used
#'
#' @return \code{data_matrix}-size data matrix with batch-effect corrected by
#'   \code{ComBat}
#'
#' @examples 
#' combat_corrected_matrix <- correct_with_ComBat(
#' example_proteome_matrix, example_sample_annotation)
#' 
#' @export
#'
correct_with_ComBat <- function(df_long, sample_annotation = NULL, 
                                feature_id_col = 'FullPeptideName',
                                measure_col = 'Intensity',
                                sample_id_col = 'FullRunName',
                                batch_col = 'MS_batch', 
                                par.prior = TRUE){
  
  data_matrix = long_to_matrix(df_long, feature_id_col = feature_id_col,
                               measure_col = measure_col, sample_id_col = sample_id_col)  
  
  #TODO: substitute with "check for sample consistency"
  sampleNames = colnames(data_matrix)
  s_a <- sample_annotation[[sample_id_col]]
  all <- union(sampleNames, s_a)
  non_matched <- all[!all %in% intersect(sampleNames, s_a)]
  if(length(non_matched)!=0){warning("Sample ID in data matrix and 
                                       sample annotation don't match. 
                                       Non-matching elements are removed for analysis")}
  
  sample_annotation = sample_annotation %>%
    filter(!!(as.name(sample_id_col)) %in% sampleNames) %>%
    arrange(match(!!(as.name(sample_id_col)), sampleNames)) %>%
    droplevels()
  
  batches = sample_annotation[[batch_col]]
  modCombat = model.matrix(~1, data = sample_annotation)
  corrected_proteome = ComBat(dat = data_matrix, batch = batches,
                                   mod = modCombat, par.prior = par.prior)
  return(corrected_proteome)
}


#' Batch correction method allows correction of 
#' continuous signal drift within batch (if required) and adjustment for
#' discrete difference across batches. 
#'
#' @inheritParams proBatch
#' @param continuous_func function to use for the fit (currently 
#' only \code{loess_regression} available); if order-associated fix is not required, should be `NULL`
#' @param discrete_func function to use for adjustment of discrete batch effects (\code{MedianCentering} or \code{ComBat})
#' @param abs_threshold the absolute threshold to filter data for curve fitting 
#' @param pct_threshold the percentage threshold to filter data for curve fitting 
#' @param ... other parameters, usually of \code{normalize_custom_fit}, and \code{fit_func}
#'
#' @return \code{data_matrix}-size data matrix with batch-effect 
#' corrected by fit and discrete functions
#' 
#' @examples 
#' batch_corrected_matrix <- correct_batch_effects(example_proteome_matrix, example_sample_annotation, 
#' continuos_func = 'loess_regression',
#' discreteFunc = 'MedianCentering', 
#' batch_col = 'MS_batch',  
#' span = 0.7,
#' abs_threshold = 5, pct_threshold = 0.20)
#' 
#' @export
#' 
#' @seealso \code{\link{fit_nonlinear}}, \code{\link{center_peptide_batch_medians}},
#' \code{\link{correct_with_ComBat}}, \code{\link{adjust_batch_trend}}
correct_batch_effects <- function(df_long, sample_annotation, 
                                  continuous_func = NULL, 
                                  discrete_func = c("MedianCentering", "ComBat"), 
                                  batch_col = 'MS_batch',  
                                  feature_id_col = 'peptide_group_label', 
                                  sample_id_col = 'FullRunName',
                                  measure_col = 'Intensity',  
                                  order_col = 'order', 
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
    median_long = center_peptide_batch_medians(df_long = df_long, 
                                               sample_annotation = sample_annotation,
                                               sample_id_col = sample_id_col,
                                               batch_col = batch_col,
                                               feature_id_col = feature_id_col,
                                               measure_col = measure_col)
    normalized_matrix = long_to_matrix(median_long, feature_id_col = feature_id_col,
                                       measure_col = measure_col, 
                                       sample_id_col = sample_id_col)
  }
  
  if(discrete_func == 'ComBat'){
    #TODO: pick up the functional programming argument borrowing
    normalized_matrix = correct_with_ComBat(df_long = df_long, 
                                            sample_annotation = sample_annotation,
                                            batch_col = batch_col, par.prior = TRUE)
  }
  
  return(normalized_matrix)
}

