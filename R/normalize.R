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

#' Median normalization of the data (per batch median)
#'
#' @name normalize
#'
#' @return
#' @export
#'
#' @examples
normalize_medians_batch <- function(df_long, sample_annotation = NULL,
                                    sample_id_col = 'FullRunName',
                                    batch_col = 'MS_batch.final',
                                    feature_id_col = 'peptide_group_label',
                                    measure_col = 'Intensity'){
  if (!(sample_id_col %in% names(df_long) & batch_col %in% names(df_long)) &
      !is.null(sample_annotation)){
    df_long = df_long %>% merge(sample_annotation)
  }
  df_normalized = df_long %>%
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
                                 batch_col = 'MS_batch',
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
  batch_table <- as.data.frame(table(sample_annotation[[batch_col]], dnn = list(batch_col)), responseName = "batch_total")
  
  
  df_normalized = df_long %>%
    filter(!is.na(UQ(as.name(measure_col)))) %>% #filter(!is.na(Intensity))
    merge(sample_annotation) %>%
    arrange_(feature_id_col, sample_order_col) %>%
    group_by_at(vars(one_of(c(feature_id_col, batch_col)))) %>% #group_by(peptide_group_label, MS_batch.final) 
    #filter(n() >3)%>%
    nest() %>%
    full_join(batch_table, by = batch_col) %>%
    mutate(fit = map2(data, batch_total, fit_func, response.var = measure_col, 
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
#' @name normalize
#' @param par.prior
#'
#' @return `data_matrix`-size data matrix with batch-effect corrected by
#'   `ComBat`
#' @export
#'
#' @examples
correct_with_ComBat <- function(data_matrix, sample_annotation, 
                                sample_id_col = 'FullRunName',
                                batch_col = 'MS_batch.final', 
                                par.prior = TRUE){

  sampleNames = colnames(data_matrix)
  sample_annotation = sample_annotation %>%
    filter(UQ(as.name(sample_id_col)) %in% sampleNames) %>%
    arrange(match(UQ(as.name(sample_id_col)), sampleNames))
  
  batches = sample_annotation[[batch_col]]
  modCombat = model.matrix(~1, data = sample_annotation)
  corrected_proteome = sva::ComBat(dat = data_matrix, batch = batches,
                              mod = modCombat, par.prior = par.prior)
  return(corrected_proteome)
}

#' Batch correction method allows correction of continuous sigal drift within batch and 
#' discrete difference across batches. 
#'
#' @name correct_batch
#' @param fitFunct function to use for the fit (currently only `loess_regression` available)
#' @param discreteFunc function to use fo discrete batch correction (`MedianCentering` or `ComBat`)
#' @param ... other parameters, usually of `normalize_custom_fit`, and `fit_func`
#'
#' @return `data_matrix`-size data matrix with batch-effect corrected by fit and discrete functions
#' @export
#'
#' @examples
correct_batch <- function(data_matrix, sample_annotation, fitFunc = 'loess_regression', 
                          discreteFunc = 'MedianCentering', batch_col = 'MS_batch',  
                          feature_id_col = 'peptide_group_label', sample_id_col = 'FullRunName',
                          measure_col = 'Intensity',  sample_order_col = 'order',...){
  
  fit_list = normalize_custom_fit(data_matrix, sample_annotation = sample_annotation,
                                  batch_col = batch_col,
                                  feature_id_col = feature_id_col,
                                  sample_id_col = sample_id_col,
                                  measure_col = measure_col,
                                  sample_order_col = sample_order_col,
                                  fit_func = fit_nonlinear,
                                  fitFunc = fitFunc)
  fit_matrix = fit_list$data_matrix
  fit_long = matrix_to_long(fit_matrix, feature_id_col = feature_id_col,
                            measure_col = measure_col, sample_id_col = sample_id_col)
  
  if(discreteFunc == 'MedianCentering'){
    median_long = normalize_medians_batch(df_long = fit_long, sample_annotation = sample_annotation,
                                          sample_id_col = sample_id_col,
                                          batch_col = batch_col,
                                          feature_id_col = feature_id_col,
                                          measure_col = measure_col)
    normalized_matrix = convert_to_matrix(median_long, feature_id_col = feature_id_col,
                                          measure_col = measure_col, sample_id_col = sample_id_col)
  }
  
  if(discreteFunc == 'ComBat'){
    batches = unique(sample_annotation[[batch_col]])
    matrix_batch = list()
    for(i in 1:length(batches)){
      samples = sample_annotation[[sample_id_col]][sample_annotation[[batch_col]] == batches[[i]]]
      matrix = fit_matrix[,samples]
      matrix_batch[[i]] = matrix[rowSums(is.na(matrix)) != ncol(matrix), ]
    }
    
    df = Reduce(function(x, y) merge(x, y, all = FALSE), 
                lapply(matrix_batch, function(y) data.table(y, keep.rownames=TRUE, key = "rn")))
    
    filtered_matrix <-as.matrix(df[,-1])
    rownames(filtered_matrix)<-unlist(df[,1])
    
    nfiltered = nrow(fit_matrix) - nrow(filtered_matrix)
    if(nfiltered > 0){
      warning(sprintf("%i rows have no measurement for one or more batches and are removed for ComBat batch correction", nfiltered))
    }
    normalized_matrix = correct_with_ComBat(filtered_matrix, sample_annotation = sample_annotation,
                                            batch_col = batch_col, par.prior = TRUE)
  }
  
  return(normalized_matrix)
}

