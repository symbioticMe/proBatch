#' Median centering of the peptides (per batch median)
#'
#' @param df_long data frame where each row is a single feature in a single
#'   sample. It minimally has a \code{sample_id_col}, a 
#'   \code{feature_id_col} and a \code{measure_col}, but 
#'   usually also an \code{m_score} (in OpenSWATH output result file)
#' @param sample_annotation data frame with sample ID, technical (e.g. MS batches) 
#'  and biological (e.g. Diet) covariates 
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if \code{df_long} is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @param batch_col column in \code{sample_annotation} that should be 
#'  used for batch comparison
#' @param feature_id_col name of the column with feature/gene/peptide/protein ID 
#' used in the long format representation df_long. In the wide formatted 
#' representation data_matrix this corresponds to the row names.
#' 
#' @return \code{df_long}-size long format data with batch-effect corrected with
#'   per-feature batch median centering in Intensity_normalized column
#'   
#' @examples 
#' median_centered_proteome <- center_peptide_batch_medians(
#' example_proteome, example_sample_annotation)
#' 
#' @export
#'
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
        mutate(median_batch = median(UQ(sym(measure_col)), na.rm = TRUE)) %>%
        ungroup() %>%
        group_by_at(vars(one_of(feature_id_col))) %>%
        mutate(median_global = median(UQ(sym(measure_col)), na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(diff = median_global - median_batch) %>%
        mutate(Intensity_normalized = UQ(sym(measure_col))+diff)
    
    return(df_normalized)
}


#' adjust batch signal trend with the custom (continuous) fit
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data frame with sample ID, technical (e.g. MS batches) 
#'  and biological (e.g. Diet) covariates 
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if \code{df_long} is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @param batch_col column in \code{sample_annotation} that should be 
#'  used for batch comparison
#' @param sample_order_col column, determining the order of sample MS run, used
#'   as covariate to fit the non-linear fit
#' @param feature_id_col name of the column with feature/gene/peptide/protein ID 
#'  used in the long format representation df_long. In the wide formatted 
#'  representation data_matrix this corresponds to the row names.
#' @param fit_func function to fit the (non)-linear trend
#' @param abs_threshold the absolute threshold to filter data for curve fitting 
#' @param pct_threshold the percentage threshold to filter data for curve fitting 
#' @param ... other parameters, usually those of the \code{fit_func}
#'
#' @return list of two items: 1) \code{data_matrix}, adjusted with continious fit; 
#' 2) fit_df, used to examine the fitting curves
#' @examples 
#' adjust_batch_trend(example_proteome_matrix, 
#' example_sample_annotation, span = 0.7, 
#' abs_threshold = 5, pct_threshold = 0.20)
#' 
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
                               abs_threshold = 5, pct_threshold = 0.20, ...){
  
  sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
  sampleNames <- colnames(data_matrix)
  s_a <- sample_annotation[[sample_id_col]]
  all <- union(sampleNames, s_a)
  non_matched <- all[!all %in% intersect(sampleNames, s_a)]
  if(length(non_matched)!=0){warning("Sample ID in data matrix 
                                     and sample annotation don't match. 
                                     Non-matching elements are removed for analysis")}
  
  sample_annotation = sample_annotation %>%
    filter(UQ(as.name(sample_id_col)) %in% sampleNames) %>%
    arrange(match(UQ(as.name(sample_id_col)), sampleNames)) %>%
    droplevels()
  
  data_matrix = as.data.frame(data_matrix)
  data_matrix[[feature_id_col]] = rownames(data_matrix)
  
  df_long = data_matrix %>%
    melt(id.vars = feature_id_col)
  names(df_long) = c(feature_id_col, sample_id_col, measure_col)
  batch_table <- as.data.frame(table(sample_annotation[[batch_col]], 
                                     dnn = list(batch_col)), 
                               responseName = "batch_total")
  sample_annotation = sample_annotation %>%
    full_join(batch_table, by = batch_col)
  
  df_normalized = df_long %>%
    filter(!is.na(UQ(as.name(measure_col)))) %>% #filter(!is.na(Intensity))
    merge(sample_annotation, by = sample_id_col) %>%
    arrange_(feature_id_col, sample_order_col) %>%
    group_by_at(vars(one_of(c(feature_id_col, batch_col, "batch_total")))) %>%  
    nest() %>%
    mutate(fit = map2(data, batch_total, fit_func, response.var = measure_col, 
                      expl.var = sample_order_col, 
                      abs_threshold = abs_threshold, 
                      pct_threshold = pct_threshold, ...)) %>%
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


#' Adjusts for discrete batch effects using ComBat
#'
#' @description Standardized input-output ComBat normalization ComBat allows users to adjust
#' for batch effects in datasets where the batch covariate is known, using
#' methodology described in Johnson et al. 2007. It uses either parametric or
#' non-parametric empirical Bayes frameworks for adjusting data for batch
#' effects. Users are returned an expression matrix that has been corrected for
#' batch effects. The input data are assumed to be cleaned and normalized before
#' batch effect removal.
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data frame with sample ID, technical (e.g. MS batches) 
#'  and biological (e.g. Diet) covariates 
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param batch_col column in \code{sample_annotation} that should be 
#'  used for batch comparison
#' @param par.prior whether parametrical or non-parametrical prior should be used
#'
#' @return \code{data_matrix}-size data matrix with batch-effect corrected by
#'   \code{ComBat}
#'
#' @examples 
#' correct_with_ComBat(example_proteome_matrix, example_sample_annotation)
#' 
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
    if(length(non_matched)!=0){warning("Sample ID in data matrix and 
                                       sample annotation don't match. 
                                       Non-matching elements are removed for analysis")}
    
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


#' Batch correction method allows correction of 
#' continuous sigal drift within batch and 
#' discrete difference across batches. 
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data frame with sample ID, technical (e.g. MS batches) 
#'  and biological (e.g. Diet) covariates 
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if \code{df_long} is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @param batch_col column in \code{sample_annotation} that should be 
#'  used for batch comparison
#' @param feature_id_col name of the column with feature/gene/peptide/protein ID 
#' used in the long format representation df_long. In the wide formatted 
#' representation data_matrix this corresponds to the row names. 
#' @param sample_order_col column, determining the order of sample MS run, used
#'   as covariate to fit the non-linear fit
#' @param fitFunc function to use for the fit (currently 
#' only \code{loess_regression} available)
#' @param discreteFunc function to use for discrete 
#' batch correction (\code{MedianCentering} or \code{ComBat})
#' @param abs_threshold the absolute threshold to filter data for curve fitting 
#' @param pct_threshold the percentage threshold to filter data for curve fitting 
#' @param ... other parameters, usually of \code{normalize_custom_fit}, and \code{fit_func}
#'
#' @return \code{data_matrix}-size data matrix with batch-effect 
#' corrected by fit and discrete functions
#' 
#' @examples 
#' correct_batch_effects(example_proteome_matrix, example_sample_annotation, 
#' discreteFunc = 'MedianCentering', 
#' batch_col = 'MS_batch',  
#' span = 0.7,
#' abs_threshold = 5, pct_threshold = 0.20)
#' 
#' @export
#'                                                 
correct_batch_effects <- function(data_matrix, sample_annotation, 
                                  fitFunc = 'loess_regression', 
                                  discreteFunc = c("MedianCentering", "ComBat"), 
                                  batch_col = 'MS_batch',  
                                  feature_id_col = 'peptide_group_label', 
                                  sample_id_col = 'FullRunName',
                                  measure_col = 'Intensity',  
                                  sample_order_col = 'order', 
                                  abs_threshold = 5, pct_threshold = 0.20, ...){
 
    discreteFunc <- match.arg(discreteFunc)
  
    sample_annotation[[batch_col]] <- as.factor(sample_annotation[[batch_col]])
    fit_list = adjust_batch_trend(data_matrix, 
                                  sample_annotation = sample_annotation,
                                  batch_col = batch_col,
                                  feature_id_col = feature_id_col,
                                  sample_id_col = sample_id_col,
                                  measure_col = measure_col,
                                  sample_order_col = sample_order_col,
                                  fit_func = fit_nonlinear,
                                  fitFunc = fitFunc, 
                                  abs_threshold = abs_threshold, 
                                  pct_threshold = pct_threshold, ...)
    fit_matrix = fit_list$data_matrix
    fit_long = matrix_to_long(fit_matrix, feature_id_col = feature_id_col,
                              measure_col = measure_col, sample_id_col = sample_id_col)
    
    if(discreteFunc == 'MedianCentering'){
        median_long = center_peptide_batch_medians(df_long = fit_long, 
                                                   sample_annotation = sample_annotation,
                                                   sample_id_col = sample_id_col,
                                                   batch_col = batch_col,
                                                   feature_id_col = feature_id_col,
                                                   measure_col = measure_col)
        normalized_matrix = long_to_matrix(median_long, feature_id_col = feature_id_col,
                                           measure_col = measure_col, 
                                           sample_id_col = sample_id_col)
    }
    
    if(discreteFunc == 'ComBat'){
        fit_matrix = long_to_matrix(fit_long, feature_id_col = feature_id_col,
                                    measure_col = measure_col, sample_id_col = sample_id_col)
        normalized_matrix = correct_with_ComBat(fit_matrix, 
                                                sample_annotation = sample_annotation,
                                                batch_col = batch_col, par.prior = TRUE)
    }
    
    return(normalized_matrix)
}

