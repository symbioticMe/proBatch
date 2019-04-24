
#' Fit a non-linear trend (currently optimized for LOESS)
#'
#' @param df_peptide data frame containing response variable e.g. 
#' samples in order and explanatory 
#'   variable e.g. measurement
#' @param batch.size the total number of samples in the batch to 
#' compute for percentage threshold 
#' @param response.var the name of the column in \code{df_peptide} (or, generally, in \code{df_long}) with 
#' the response variable, usually \code{measure_col}
#' @param expl.var the name of the column in \code{df_peptide} with the 
#' explanatory variable, in most cases order
#' @param noFitRequants (logical) whether to fit requanted values
#' @param fitFunc function to use for the fit (\code[loess_regression})
#' @param optimize_span logical, whether to specify span or optimize it
#' @param qual_col column, indicating whether the "quality" variable is specified
#' @param qual_value 
#' @param abs_threshold the absolute threshold to filter 
#' data for curve fitting 
#' @param pct_threshold the percentage threshold to filter 
#' data for curve fitting 
#' @param ... additional paramters to be passed to the fitting function
#'
#' @return vector of fitted response values
#' 
#' @keywords internal
#' 
fit_nonlinear <- function(df_peptide, batch.size = NULL, response.var = 'y', expl.var = 'x',
                          feature.id = NULL, batch.id = NULL,
                          fitFunc = 'loess_regression',
                          optimize_span = FALSE, 
                          noFitRequants = FALSE, qual_col = 'm_score', qual_value = 2,
                          abs_threshold = 5, pct_threshold = 0.20, ...){

    df_peptide <- df_peptide[sort.list(df_peptide[[expl.var]]),]
    x_all = df_peptide[[expl.var]]
    y = df_peptide[[response.var]]
    
    if(noFitRequants){
      if(!is.null(qual_col) && (qual_col %in% names(df_peptide))){
        warning('imputed value column is in the data, fitting curve only to measured, non-imputed values')
        imputed_values <- df_peptide[[qual_col]] == qual_value
        x_to_fit = x_all[!imputed_values]
        y = y[!imputed_values]
      } else {
        stop('imputed values are specified not to be used for curve fitting, however, no flag for imputed values is specified')
      }
      
    } else {
      if(!is.null(qual_col) && (qual_col %in% names(df_peptide))){
        warning('requant column is in the data, are you sure you want to fit non-linear curve to these values, too?')
      }
      x_to_fit = x_all
    } 
    
    #checking if there is a reasonable number of values to fit any sensible curve
    if (is.null(batch.size)){
      warning("Batch size is not specified, assuming number of entries for this 
              batch and feature is the total batch size. Maybe erroneous,
              if missing values are not NAs, but removed from data frame completely")
      batch.size = nrow(df_peptide)
    }
    
    pct_threshold = batch.size*pct_threshold
    if(length(x_to_fit) >= abs_threshold & length(x_to_fit) >= pct_threshold){
      #fitting the curve
      #TODO: re-write in the functional programming paradigm (e.g. arguments - function, x_all, y, x_to_fit)
      if(fitFunc == 'loess_regression'){
        if(!optimize_span){
          fit_res = loess_regression(x_all, y, x_to_fit,...)
        } else {
          fit_res = loess_regression_opt(x_to_fit, y, x_all, ...)
        }
      } else{
        stop("Only loess regression fitting is available for current version")
      }
    }else{
      warning(sprintf("Curve fitting didn't have enough points to fit for the feature %s
                      in the batch %s, leaving the original value", feature.id, batch.id))
      fit_res = y
    }
  return(fit_res)
}

loess_regression <- function(x_all, y, x_to_fit, ...){
    fit = loess(y ~ x_to_fit, surface = 'direct', ...)
    pred <- predict(fit, newdata = x_all)
    return(pred)
}

loess_regression_opt <- function(x_to_fit, y, x_all, ...) {
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
  loo.values.bw <- sapply(1:n, loo.reg.value.bw, x, y, bw, kernel)
  ## Calculate and return MSE
  return(mean((y - loo.values.bw)^2))
}

optimise_bw <- function(x, y, kernel = 'normal',
                        bws = c(0.01, 0.5, 1, 1.5, 2, 5, 10)){
  cv.nw_mult = sapply(bws, function(bw) loocv.nw(x, y, bw, kernel))
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
  for (j in 1:n)
      S.nw[, j] <- reg.fcn.nw(x, Id[,j],x, bw = bw)
    df.nw <- sum(diag(S.nw))
    return(df.nw)
  }