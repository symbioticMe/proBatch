
#' Fit a non-linear trend (currently optimized for LOESS)
#'
#' @inheritParams proBatch
#' @param df_feature_batch data frame containing response variable e.g. 
#' samples in order and explanatory 
#'   variable e.g. measurement for a specific feature (peptide) in a specific batch
#' @param batch_size the total number of samples in the batch to 
#' compute for percentage threshold 
#' @param feature_id the name of the feature, required for warnings
#' @param batch_id  the name of the batch, required for warnings
#' @param fit_func function to use for the fit, e.g. \code{loess_regression}
#' @param optimize_span logical, whether to specify span or optimize it 
#' (specific entirely for LOESS regression)
#' @param no_fit_imputed (logical) whether to fit the imputed (requant) values
#' @param abs_threshold the absolute threshold to filter 
#' data for curve fitting 
#' @param pct_threshold the percentage threshold to filter 
#' data for curve fitting
#' 
#' @param ... additional parameters to be passed to the fitting function
#'
#' @return vector of fitted response values
#' 
#' @examples 
#' test_peptide = example_proteome$peptide_group_label[1] 
#' df_selected = example_proteome[example_proteome$peptide_group_label == test_peptide,]
#' batch_selected_df = example_sample_annotation[example_sample_annotation$MS_batch == 'Batch_1',]
#' df_for_test = merge(df_selected, batch_selected_df, by = 'FullRunName')
#' fit_values = fit_nonlinear(df_for_test)
#' 
#' @export
#' 
fit_nonlinear <- function(df_feature_batch, batch_size = NULL, 
                          measure_col = 'Intensity', order_col = 'order',
                          feature_id = NULL, batch_id = NULL,
                          fit_func = 'loess_regression',
                          optimize_span = FALSE, 
                          no_fit_imputed = FALSE, qual_col = 'm_score', qual_value = 2,
                          abs_threshold = 5, pct_threshold = 0.20, ...){
  #df_feature_batch <- df_feature_batch[sort.list(df_feature_batch[[order_col]]),]
  x_all = df_feature_batch[[order_col]]
  y = df_feature_batch[[measure_col]]
    
  if(no_fit_imputed){
    if(!is.null(qual_col) && (qual_col %in% names(df_feature_batch))){
      warning('imputed value column is in the data, fitting curve only to measured, non-imputed values')
      imputed_values <- df_feature_batch[[qual_col]] == qual_value
      x_to_fit = x_all[!imputed_values]
      y[imputed_values] = NA
    } else {
      stop('imputed values are specified not to be used for curve fitting, however, 
           no flag for imputed values is specified')
      }
    } else {
      if(!is.null(qual_col) && (qual_col %in% names(df_feature_batch))){
        warning('imputed value (requant) column is in the data, are you sure you want to fit non-linear curve to these values, too?')
      }
      x_to_fit = x_all
    }
    
    #checking if there is a reasonable number of values to fit any sensible curve
    if (is.null(batch_size)){
      warning("Batch size is not specified, assuming number of entries for this 
              batch and feature is the total batch size. Maybe erroneous,
              if missing values are not NAs, but removed from data frame completely")
      batch_size = nrow(df_feature_batch)
    }
    
    pct_threshold = batch_size*pct_threshold
    if(length(x_to_fit) >= abs_threshold & length(x_to_fit) >= pct_threshold){
      #fitting the curve
      #TODO: re-write in the functional programming paradigm (e.g. arguments - function, x_all, y, x_to_fit)
      if(fit_func == 'loess_regression'){
        if(!optimize_span){
          fit_res = loess_regression(x_to_fit, y, x_all, 
                                     feature_id = feature_id, batch_id = batch_id, ...)
        } else {
          fit_res = loess_regression_opt(x_to_fit, y, x_all, 
                                         feature_id = feature_id, batch_id = batch_id, ...)
        }
      } else{
        stop("Only loess regression fitting is available for current version")
      }
    }else{
      warning(sprintf("Curve fitting didn't have enough points to fit for the feature %s
                      in the batch %s, leaving the original value", feature_id, batch_id))
      fit_res = y
    }
  fit_res[is.na(y)] = NA
  return(fit_res)
}

loess_regression <- function(x_to_fit, y, x_all,   
                             feature_id = NULL, batch_id = NULL, ...){
  out <- tryCatch({
    fit = loess(y ~ x_to_fit, surface = 'direct', ...)
    pred <- predict(fit, newdata = x_all)
    pred
  },
    warning=function(cond) {
      message(sprintf("Feature %s in batch %s caused a warning:", feature_id, batch_id))
      message("Here's the original warning message:")
      message(cond)
      # Choose a return value in case of warning
      message(sprintf('class of prediction is %s'), class(pred))
      message(sprintf('length of prediction is %s'), length(pred))
      y
    }
  )
  return(out)
}

loess_regression_opt <- function(x_to_fit, y, x_all, 
                                 feature_id = NULL, batch_id = NULL, ...) {
  bw = optimise_bw(x_to_fit, y, ...)
  degr_freedom = optimise_df(x_to_fit, bw)
  fit = loess(y ~ x_all, enp.target = degr_freedom, surface = 'direct', ...)
  res = predict(fit, newdata = x_to_fit)
  return(res)
}

loocv.nw <- function(x, y, bw = 1.5, kernel = "normal"){
  ## Help function to calculate leave-one-out regression values
  loo.reg.value.bw <- function(i, x, y, bw, kernel = kernel){
    return(ksmooth(x[-i], y[-i],
                   x.points = x[i],
                   kernel = kernel, bandwidth = bw)$y)
  }
    
    
  ## Calculate LOO regression values using the help function above
  
  n <- max(length(x), length(y))
  loo.values.bw <- vapply(seq_len(n), FUN = loo.reg.value.bw, 
                          FUN.VALUE = numeric(1),
                          x, y, bw, kernel)
  ## Calculate and return MSE
  return(mean((y - loo.values.bw)^2))
}

optimise_bw <- function(x, y, kernel = 'normal',
                        bws = c(0.01, 0.5, 1, 1.5, 2, 5, 10)){
  cv.nw_mult = vapply(bws, FUN = function(bw) loocv.nw(x, y, bw, kernel), 
                      FUN.VALUE = numeric(1))
  bw_best = bws[which.min(cv.nw_mult)]
  return(bw_best)
}

reg.fcn.nw <- function(reg.x, reg.y, x, bw = 1.5)
  ksmooth(reg.x, reg.y, x.points = x, kernel = "normal", bandwidth = bw)$y

optimise_df <- function(x, bw){
  #find the best bandwidth
  
  n <- length(x)
  Id <- diag(n)
  S.nw <- matrix(0, n, n)
  for (j in seq_len(n))
      S.nw[, j] <- reg.fcn.nw(x, Id[,j],x, bw = bw)
    df.nw <- sum(diag(S.nw))
    return(df.nw)
  }
