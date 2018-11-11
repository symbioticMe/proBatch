#' Fit a non-linear trend
#'
#' @param dataDF
#' @param batch.size the total number of samples in the batch to compute for percentage threshold 
#' @param response.var the name of the column in dataDF with the response variable
#' @param expl.var the name of the column in dataDF with the explanatory variable
#' @param noFitRequants (logical) whether to fit requanted values
#' @param fitFunc function to use for the fit (`loess_regression`)
#' @param with_df
#' @param loess.span the parameter α which controls the degree of smoothing for loess
#' @param abs.threshold the absolute threshold to filter data for curve fitting 
#' @param pct.threshold the percentage threshold to filter data for curve fitting 
#' @param ... additional paramters to be passed to the fitting function
#'
#' @return vector of fitted response values
#'
#' @export
#' @keywords internal

# TODO: Document dataDF and with_df
fit_nonlinear <- function(dataDF, batch.size, response.var = 'y', expl.var = 'x',
                          noFitRequants = F, fitFunc = 'loess_regression',
                          with_df = F, loess.span = 0.75, abs.threshold = 5, pct.threshold = 0.20, ...){
  
  dataDF <- dataDF[sort.list(dataDF[[expl.var]]),]
  x_to_fit = dataDF[[expl.var]]
  y = dataDF[[response.var]]
  pct.threshold = batch.size*pct.threshold
  if(fitFunc == "loess_regression"){
    if(length(x_to_fit)*loess.span >= abs.threshold & length(x_to_fit)*loess.span >= pct.threshold){
      x_all = x_to_fit
      if(noFitRequants){
        x_all[dataDF$requant] = NA
      }
      if(!with_df){
        fit_res = switch(fitFunc,
                         loess_regression = loess_regression(x_all, y, x_to_fit, span = loess.span,...)
        )
      } else {
        bw = optimise_bw(dataDF, response.var = response.var, expl.var = expl.var)
        df = optimise_df(dataDF, bw, response.var = response.var, expl.var = expl.var)
        fit_res = switch(fitFunc,
                         loess_regression = loess_regression_opt(x_all, y, x_to_fit, df, span = loess.span,...))
      }
    }else{
      fit_res = rep(NA, length(x_to_fit))
    }
  }else{
   stop("Only loess regression fitting is available for current version")
  }
  return(fit_res)
}

loess_regression <- function(x_all, y, x_to_fit, ...){
  fit = loess(y ~ x_to_fit, surface = 'direct', ...)
  pred <- predict(fit, newdata = x_to_fit)
  return(pred)
}

loess_regression_opt <- function(x_all, y, x_to_fit, df, ...){
  fit = loess(y ~ x_all, enp.target = df, surface = 'direct', ...)
  res = predict(fit, newdata = x_to_fit)
  return(res)
}
