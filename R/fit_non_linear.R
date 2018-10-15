# TODO: Make sure all options work

#' Fit a non-linear trend
#'
#' Fit a non-linear trend
#'
#' @param dataDF
#' @param response.var the name of the column in dataDF with the response variable
#' @param expl.var the name of the column in dataDF with the explanatory variable
#' @param noFitRequants (logical) whether to fit requanted values
#' @param fitFunc function to use for the fit (`kernel_smooth`, `smooth_spline`, or `loess_regression`)
#' @param with_df
#' @param loess.span the parameter α which controls the degree of smoothing for loess
#' @param ... additional paramters to be passed to the fitting function
#'
#' @return vector of fitted response values
#'
#' @export
#' @keywords internal

# TODO: Document dataDF and with_df
fit_nonlinear <- function(dataDF, response.var = 'y', expl.var = 'x',
                          noFitRequants = F, fitFunc = 'kernel_smooth',
                          with_df = F, loess.span = 0.75, ...){
  dataDF <- dataDF[sort.list(dataDF[[expl.var]]),]
  x_to_fit = dataDF[[expl.var]]
  y = dataDF[[response.var]]
  if(fitFunc == "loess_regression"){
    if(length(x_to_fit)*loess.span >= 3){
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
      x_all = x_to_fit
      if(noFitRequants){
        x_all[dataDF$requant] = NA
      }
      if(!with_df){
        fit_res = switch(fitFunc,
                         kernel_smooth = kernel_smooth(x_all, y, x_to_fit, ...),
                         smooth_spline = smooth_spline(x_all, y, x_to_fit, ...))
      } else {
        bw = optimise_bw(dataDF, response.var = response.var, expl.var = expl.var)
        df = optimise_df(dataDF, bw, response.var = response.var, expl.var = expl.var)
        fit_res = switch(fitFunc,
                         kernel_smooth = kernel_smooth_opt(x_all, y, x_to_fit, bw, ...),
                         smooth_spline = smooth_spline_opt(x_all, y, x_to_fit, df, ...))
      }
    }
  return(fit_res)
}




loocv.nw <- function(reg.data, bw = 1.5, kernel = "normal",
                     response.var = 'y', expl.var = 'x', ...){
  ## Help function to calculate leave-one-out regression values
  loo.reg.value.bw <- function(i, reg.data, bw, kernel = kernel,
                               response.var = 'y', expl.var = 'x')
    return(ksmooth(reg.data[[expl.var]][-i], reg.data[[response.var]][-i],
                   x.points = reg.data[[expl.var]][i],
                   kernel = kernel, bandwidth = bw)$y)

  ## Calculate LOO regression values using the help function above
  n <- nrow(reg.data)
  loo.values.bw <- sapply(1:n, loo.reg.value.bw, reg.data, bw, kernel,
                          response.var, expl.var)

  ## Calculate and return MSE
  return(mean((reg.data[[response.var]] - loo.values.bw)^2))
}

optimise_bw <- function(dataDF, response.var = 'y', expl.var = 'x',
                        bws = c(0.01, 0.5, 1, 1.5, 2, 5, 10)){
  cv.nw_mult = sapply(bws, function(bw) loocv.nw(dataDF, bw,
                                                 response.var = response.var,
                                                 expl.var = expl.var))
  bw_best = bws[which.min(cv.nw_mult)]
  return(bw_best)
}

reg.fcn.nw <- function(reg.x, reg.y, x, bw = 1.5)
  ksmooth(reg.x, reg.y, x.points = x, kernel = "normal", bandwidth = bw)$y

optimise_df <- function(dataDF, bw, response.var = 'y', expl.var = 'x'){
  #find the best bandwidth

  n <- nrow(dataDF)
  Id <- diag(n)
  S.nw <- matrix(0, n, n)
  for (j in 1:n)
    S.nw[, j] <- reg.fcn.nw(dataDF[[expl.var]], Id[,j], dataDF[[expl.var]], bw = bw)
  df.nw <- sum(diag(S.nw))
  return(df.nw)
}

kernel_smooth <- function(x_all, y, x_to_fit, bw = 1.5, ...){
  res = ksmooth(x_to_fit, y, x.point = x_all,
          kernel = "normal", bandwidth = bw)$y
  return(res)
}

smooth_spline <- function(x_all, y, x_to_fit, ...){
  #allow to define nknots among other parameters or spar or whatever

  nknots = length(x_all)/2

  fit = smooth.spline(x_to_fit, y, nknots = nknots)
  res = predict(fit, x_all)$y
  return(res)
}

loess_regression <- function(x_all, y, x_to_fit, ...){
  fit = loess(y ~ x_to_fit, surface = 'direct', ...)
  pred <- predict(fit, newdata = x_to_fit)
  return(pred)
}

kernel_smooth_opt <- function(x_all, y, x_to_fit, bw, kernel = "normal", ...){
  temp_x_all = x_all[!is.na(x_all)]
  y_for_model = y[!is.na(x_all)]
  res = ksmooth(temp_x_all, y_for_model, x.points = x_to_fit,
                kernel = kernel, bandwidth = bw)$y
  return(res)
}

smooth_spline_opt <- function(x_all, y, x_to_fit, df,...){
  temp_x_all = x_all[!is.na(x_all)]
  y_for_model = y[!is.na(x_all)]
  fit = smooth.spline(x = temp_x_all, y = y_for_model, df = df)
  res = predict(fit, x = x_to_fit)$y
  return(res)
}

loess_regression_opt <- function(x_all, y, x_to_fit, df, ...){
  fit = loess(y ~ x_all, enp.target = df, surface = 'direct', ...)
  res = predict(fit, newdata = x_to_fit)
  return(res)
}
