#' Fit a non-linear trend (currently optimized for LOESS)
#'
#' @inheritParams proBatch
#' @param df_feature_batch data frame containing response variable e.g.
#' samples in order and explanatory variable e.g. measurement for a
#' specific feature (peptide) in a specific batch
#' @param feature_id the name of the feature, required for warnings
#' @param batch_id  the name of the batch, required for warnings
#' @param fit_func function to use for the fit, e.g. \code{loess_regression}
#' @param optimize_span logical, whether to specify span or optimize it
#' (specific entirely for LOESS regression)
#' @param no_fit_imputed (logical) whether to fit the imputed (requant) values
#' @param min_measurements the absolute threshold to filter
#'
#' @param ... additional parameters to be passed to the fitting function
#'
#' @return vector of fitted response values
#'
#' @examples
#' # Load necessary datasets
#' data(list = c("example_proteome", "example_sample_annotation"), package = "proBatch")
#'
#' test_peptide <- example_proteome$peptide_group_label[1]
#' selected_peptide <- example_proteome$peptide_group_label == test_peptide
#' df_selected <- example_proteome[selected_peptide, ]
#' selected_batch <- example_sample_annotation$MS_batch == "Batch_1"
#' batch_selected_df <- example_sample_annotation[selected_batch, ]
#' df_for_test <- merge(df_selected, batch_selected_df, by = "FullRunName")
#' fit_values <- fit_nonlinear(df_for_test)
#'
#' # for the case where are two many missing values, no curve is fit
#' selected_batch <- example_sample_annotation$MS_batch == "Batch_2"
#' batch_selected_df <- example_sample_annotation[selected_batch, ]
#' df_for_test <- merge(df_selected, batch_selected_df, by = "FullRunName")
#' fit_values <- fit_nonlinear(df_for_test)
#' missing_values <- df_for_test[["m_score"]] == 2
#' all(fit_values[!is.na(fit_values)] == df_for_test[["Intensity"]][!missing_values])
#'
#' @export
#'
fit_nonlinear <- function(df_feature_batch,
                          measure_col = "Intensity", order_col = "order",
                          feature_id = NULL, batch_id = NULL,
                          fit_func = "loess_regression",
                          optimize_span = FALSE,
                          no_fit_imputed = TRUE, qual_col = "m_score",
                          qual_value = 2,
                          min_measurements = 8, ...) {
  x_all <- df_feature_batch[[order_col]]
  y_all <- df_feature_batch[[measure_col]]

  if (no_fit_imputed) {
    if (!is.null(qual_col) && (qual_col %in% names(df_feature_batch))) {
      warning("imputed value column is in the data, fitting curve only to
              measured, non-imputed values")
      imputed_values <- df_feature_batch[[qual_col]] == qual_value
      df_feature_batch[[measure_col]][imputed_values] <- NA
      missing_values <- imputed_values
    } else {
      stop("imputed values are specified not to be used for curve fitting,
           however, no flag for imputed values is specified")
    }
  } else {
    if (!is.null(qual_col) && (qual_col %in% names(df_feature_batch))) {
      warning("imputed value (requant) column is in the data, are you sure you
                want to fit non-linear curve to these values, too?")
    }
    missing_values <- is.na(y_all)
  }

  y <- y_all[!missing_values]
  x_to_fit <- x_all[!missing_values]

  max_consec_meas <- rle_func(df_feature_batch)
  if (max_consec_meas >= min_measurements) {
    # fitting the curve
    # TODO: re-write in the functional programming paradigm (e.g. arguments -
    #       function, x_all, y, x_to_fit)
    if (fit_func == "loess_regression") {
      if (!optimize_span) {
        fit_res <- loess_regression(x_to_fit, y, x_all, y_all,
          feature_id = feature_id,
          batch_id = batch_id, ...
        )
      } else {
        fit_res <- loess_regression_opt(x_to_fit, y, x_all, y_all,
          feature_id = feature_id,
          batch_id = batch_id, ...
        )
      }
    } else {
      stop("Only loess regression fitting is available for current version")
    }
  } else {
    warning(sprintf(
      "Curve fitting didn't have enough points to fit for the
                       feature %s in the batch %s, leaving the original value",
      feature_id, batch_id
    ))
    fit_res <- rep(NA, length(y_all))
  }

  return(fit_res)
}

loess_regression <- function(x_to_fit, y, x_all, y_all,
                             feature_id = NULL, batch_id = NULL, ...) {
  out <- tryCatch(
    {
      fit <- loess(y ~ x_to_fit, surface = "direct", ...)
      pred <- predict(fit, newdata = x_all)
      pred
    },
    warning = function(cond) {
      message(sprintf(
        "Feature %s in batch %s could not be fit with LOESS:",
        feature_id, batch_id
      ))
      message(cond)
      rep(NA, length(y_all))
    }
  )
  return(out)
}

loess_regression_opt <- function(x_to_fit, y, x_all, y_all,
                                 feature_id = NULL, batch_id = NULL, ...) {
  out <- tryCatch(
    {
      bw <- optimise_bw(x_to_fit, y, ...)
      degr_freedom <- optimise_df(x_to_fit, bw)
      fit <- loess(y ~ x_all, enp.target = degr_freedom, surface = "direct", ...)
      pred <- predict(fit, newdata = x_to_fit)
      return(pred)
    },
    warning = function(cond) {
      message(sprintf(
        "Feature %s in batch %s could not be fit with optimised LOESS:",
        feature_id, batch_id
      ))
      message(cond)
      y_all
    }
  )
  return(out)
}

loocv.nw <- function(x, y, bw = 1.5, kernel = "normal") {
  ## Help function to calculate leave-one-out regression values
  loo.reg.value.bw <- function(i, x, y, bw, kernel = kernel) {
    return(ksmooth(x[-i], y[-i],
      x.points = x[i],
      kernel = kernel, bandwidth = bw
    )$y)
  }


  ## Calculate LOO regression values using the help function above

  n <- max(length(x), length(y))
  loo.values.bw <- vapply(seq_len(n),
    FUN = loo.reg.value.bw,
    FUN.VALUE = numeric(1),
    x, y, bw, kernel
  )
  ## Calculate and return MSE
  return(mean((y - loo.values.bw)^2))
}

optimise_bw <- function(x, y, kernel = "normal",
                        bws = c(0.01, 0.5, 1, 1.5, 2, 5, 10)) {
  cv.nw_mult <- vapply(bws,
    FUN = function(bw) loocv.nw(x, y, bw, kernel),
    FUN.VALUE = numeric(1)
  )
  bw_best <- bws[which.min(cv.nw_mult)]
  return(bw_best)
}

reg.fcn.nw <- function(reg.x, reg.y, x, bw = 1.5) {
  ksmooth(reg.x, reg.y, x.points = x, kernel = "normal", bandwidth = bw)$y
}

optimise_df <- function(x, bw) {
  # find the best bandwidth

  n <- length(x)
  Id <- diag(n)
  S.nw <- matrix(0, n, n)
  for (j in seq_len(n)) {
    S.nw[, j] <- reg.fcn.nw(x, Id[, j], x, bw = bw)
  }
  df.nw <- sum(diag(S.nw))
  return(df.nw)
}

rle_func <- function(df) {
  rle_res <- rle(is.na(df$Intensity[order(df$order)]))
  max_measured <- max(rle_res$lengths[!rle_res$values])
  return(max_measured)
}
