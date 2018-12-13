#' Correction of batch effects in the data
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param df_long data frame where each row is a single feature in a single
#'   sample. It minimally has a \code{sample_id_col}, a \code{feature_id_col} and a
#'   \code{measure_col}, but usually also an \code{m_score} (in OpenSWATH output result
#'   file)
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if `df_long` is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @name correct_batch
NULL
#> NULL

#' Median centering of the peptides (per batch median)
#'
#' @name correct_batch
#'
#' @return `data_matrix`-size data matrix with batch-effect corrected with
#'   per-feature batch median centering
#' @export
#'
center_peptide_batch_medians <- function(df_long, sample_annotation = NULL,
                                  sample_id_col = 'FullRunName',
                                  batch_col = 'MS_batch',
                                  feature_id_col = 'peptide_group_label',
                                  measure_col = 'Intensity'){
  
  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(df_long[[sample_id_col]])) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  
  if (!(sample_id_col %in% names(df_long) & batch_col %in% names(df_long)) &
      !is.null(sample_annotation)){
    df_long = df_long %>% merge(sample_annotation, by = sample_id_col)
  }
  df_normalized = df_long %>%
    group_by_at(vars(one_of(batch_col, feature_id_col))) %>%
    mutate(median_batch = median(UQ(sym(measure_col)), na.rm = T)) %>%
    ungroup() %>%
    group_by_at(vars(one_of(feature_id_col))) %>%
    mutate(median_global = median(UQ(sym(measure_col)), na.rm = T)) %>%
    ungroup() %>%
    mutate(diff = median_global - median_batch) %>%
    mutate(Intensity_normalized = UQ(sym(measure_col))+diff)
  
  return(df_normalized)
}


#' adjust batch signal trend with the custom (continuous) fit
#'
#' @name correct_batch
#' @param sample_order_col column, determining the order of sample MS run, used
#'   as covariate to fit the non-linear fit
#' @param fit_func function to fit the (non)-linear trend
#' @param return_long whether the result should be the "long" data frame (as
#'   `df_long`) or "wide" (as `data_matrix`)
#' @param ... other parameters, usually those of the `fit_func`
#'
#' @return list of two items: 1) `data_matrix`, adjusted with continious fit; 
#' 2) fit_df, used to examine the fitting curves
#' @export
#'
#' @seealso \code{\link{fit_nonlinear}}
adjust_batch_trend <- function(data_matrix, sample_annotation,
                                 batch_col = 'MS_batch',
                                 feature_id_col = 'peptide_group_label',
                                 sample_id_col = 'FullRunName',
                                 measure_col = 'Intensity',
                                 sample_order_col = 'order',
                                 fit_func = fit_nonlinear, 
                                 abs.threshold = 5, pct.threshold = 0.20, ...){
  
  sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
  sampleNames <- colnames(data_matrix)
  s_a <- sample_annotation[[sample_id_col]]
  all <- union(sampleNames, s_a)
  non_matched <- all[!all %in% intersect(sampleNames, s_a)]
  if(length(non_matched)!=0){warning("Sample ID in data matrix and sample annotation don't match. Non-matching elements are removed for analysis")}
  
  sample_annotation = sample_annotation %>%
    filter(UQ(as.name(sample_id_col)) %in% sampleNames) %>%
    arrange(match(UQ(as.name(sample_id_col)), sampleNames)) %>%
    droplevels()
  
  data_matrix = as.data.frame(data_matrix)
  data_matrix[[feature_id_col]] = rownames(data_matrix)
  
  df_long = data_matrix %>%
    melt(id.vars = feature_id_col)
  names(df_long) = c(feature_id_col, sample_id_col, measure_col)
  batch_table <- as.data.frame(table(sample_annotation[[batch_col]], dnn = list(batch_col)), responseName = "batch_total")
  sample_annotation = sample_annotation %>%
    full_join(batch_table, by = batch_col)
  
  df_normalized = df_long %>%
    filter(!is.na(UQ(as.name(measure_col)))) %>% #filter(!is.na(Intensity))
    merge(sample_annotation, by = sample_id_col) %>%
    arrange_(feature_id_col, sample_order_col) %>%
    group_by_at(vars(one_of(c(feature_id_col, batch_col, "batch_total")))) %>% #group_by(peptide_group_label, MS_batch.final, tota_batch) 
    nest() %>%
    mutate(fit = map2(data, batch_total, fit_func, response.var = measure_col, 
                      expl.var = sample_order_col, 
                      abs.threshold = abs.threshold, pct.threshold = pct.threshold, ...)) %>%
    unnest() %>%
    group_by_at(vars(one_of(c(feature_id_col, batch_col)))) %>%
    mutate(mean_fit = mean(fit)) %>%
    mutate(diff = mean_fit - fit) %>%
    mutate_(Intensity_normalized = interp(~`+`(x, y),
                                          x = as.name('diff'),
                                          y = as.name(measure_col)))
  
  fit_df = df_normalized %>% dplyr::select(one_of(c('fit', feature_id_col,
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
#' @name correct_batch
#' @param par.prior whether parametrical or non-parametrical prior should be used
#'
#' @return `data_matrix`-size data matrix with batch-effect corrected by
#'   `ComBat`
#' @export
#'
correct_with_ComBat <- function(data_matrix, sample_annotation, 
                                sample_id_col = 'FullRunName',
                                batch_col = 'MS_batch', 
                                par.prior = TRUE){
  
  sampleNames = colnames(data_matrix)
  s_a <- sample_annotation[[sample_id_col]]
  all <- union(sampleNames, s_a)
  non_matched <- all[!all %in% intersect(sampleNames, s_a)]
  if(length(non_matched)!=0){warning("Sample ID in data matrix and sample annotation don't match. Non-matching elements are removed for analysis")}
  
  sample_annotation = sample_annotation %>%
    filter(UQ(as.name(sample_id_col)) %in% sampleNames) %>%
    arrange(match(UQ(as.name(sample_id_col)), sampleNames)) %>%
    droplevels()
  
  batches = sample_annotation[[batch_col]]
  modCombat = model.matrix(~1, data = sample_annotation)
  corrected_proteome = sva::ComBat(dat = data_matrix, batch = batches,
                                   mod = modCombat, par.prior = par.prior)
  return(corrected_proteome)
}


#' Batch correction method allows correction of continuous sigal drift within batch and 
#' discrete difference across batches. 
#'
#' @name correct_batch_trend
#' @param fitFunct function to use for the fit (currently only `loess_regression` available)
#' @param discreteFunc function to use for discrete batch correction (`MedianCentering` or `ComBat`)
#' @param ... other parameters, usually of `normalize_custom_fit`, and `fit_func`
#'
#' @return `data_matrix`-size data matrix with batch-effect corrected by fit and discrete functions
#' @export
#'
#' @examples \donotrun{
#' batch_corrected_matrix <- correct_batch_effects(data_matrix = quantile_normalized_matrix, 
#'                                                 example_sample_annotation, discreteFunc = 'ComBat',
#'                                                 batch_col = "MS_batch", span = 0.8,
#'                                                 abs.threshold = 5, pct.threshold = 0.20)
#' }
correct_batch_effects <- function(data_matrix, sample_annotation, fitFunc = 'loess_regression', 
                                discreteFunc = 'MedianCentering', batch_col = 'MS_batch',  
                                feature_id_col = 'peptide_group_label', sample_id_col = 'FullRunName',
                                measure_col = 'Intensity',  sample_order_col = 'order', 
                                loess.span = 0.75, abs.threshold = 5, pct.threshold = 0.20, ...){
  sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
  fit_list = adjust_batch_trend(data_matrix, sample_annotation = sample_annotation,
                                  batch_col = batch_col,
                                  feature_id_col = feature_id_col,
                                  sample_id_col = sample_id_col,
                                  measure_col = measure_col,
                                  sample_order_col = sample_order_col,
                                  fit_func = fit_nonlinear,
                                  fitFunc = fitFunc, 
                                  abs.threshold = abs.threshold, 
                                  pct.threshold = pct.threshold, ...)
  fit_matrix = fit_list$data_matrix
  fit_long = matrix_to_long(fit_matrix, feature_id_col = feature_id_col,
                            measure_col = measure_col, sample_id_col = sample_id_col)
  
  if(discreteFunc == 'MedianCentering'){
    median_long = center_peptide_batch_medians(df_long = fit_long, sample_annotation = sample_annotation,
                                          sample_id_col = sample_id_col,
                                          batch_col = batch_col,
                                          feature_id_col = feature_id_col,
                                          measure_col = measure_col)
    normalized_matrix = long_to_matrix(median_long, feature_id_col = feature_id_col,
                                       measure_col = measure_col, sample_id_col = sample_id_col)
  }
  
  if(discreteFunc == 'ComBat'){
    filtered_long = remove_peptides_with_missing_batch(fit_long, sample_annotation,
                                                       batch_col = batch_col,
                                                       feature_id_col = feature_id_col,
                                                       sample_id_col = sample_id_col)
    filtered_matrix = long_to_matrix(filtered_long, feature_id_col = feature_id_col,
                                     measure_col = measure_col, sample_id_col = sample_id_col)
    
    nfiltered = nrow(fit_matrix) - nrow(filtered_matrix)
    if(nfiltered > 0){
      warning(sprintf("%i rows have no measurement for one or more batches and are removed for ComBat batch correction", nfiltered))
    }
    normalized_matrix = correct_with_ComBat(filtered_matrix, sample_annotation = sample_annotation,
                                            batch_col = batch_col, par.prior = TRUE)
  }
  
  return(normalized_matrix)
}
