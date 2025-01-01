#' Calculate CV distribution for each feature
#'
#' @inheritParams proBatch
#' @inheritParams transform_raw_data
#' @param biospecimen_id_col column in \code{sample_annotation}
#' that defines a unique bio ID, which is usually a
#' combination of conditions or groups.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column
#' @param unlog (logical) whether to reverse log transformation of the original data
#'
#' @return data frame with Total CV for each feature & (optionally) per-batch CV
#' @export
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome"), package = "proBatch")
#' CV_df <- calculate_feature_CV(example_proteome,
#'   sample_annotation = example_sample_annotation,
#'   measure_col = "Intensity",
#'   batch_col = "MS_batch"
#' )
calculate_feature_CV <- function(df_long, sample_annotation = NULL,
                                 feature_id_col = "peptide_group_label",
                                 sample_id_col = "FullRunName",
                                 measure_col = "Intensity", batch_col = NULL,
                                 biospecimen_id_col = NULL,
                                 unlog = TRUE, log_base = 2, offset = 1) {
  df_long <- check_sample_consistency(sample_annotation, sample_id_col, df_long)

  if (is.null(biospecimen_id_col)) {
    warning("considering all samples as replicates!")
    df_long$biospecimen_id_col <- "replication"
    biospecimen_id_col <- "biospecimen_id_col"
  } else {
    if (!(biospecimen_id_col %in% names(df_long))) {
      stop("biospecimen ID, indicating replicates, is not in the data (df_long or sample_annotation)")
    }
  }

  if (unlog) {
    warning("reversing log-transformation for CV calculation!")
    df_long <- unlog_df(df_long, log_base = log_base, offset = offset, measure_col = measure_col)
  }


  if (!is.null(batch_col)) {
    df_long <- df_long %>%
      group_by(!!!syms(c(feature_id_col, batch_col, biospecimen_id_col))) %>%
      mutate(n = sum(!is.na(!!sym(measure_col))))
  } else {
    df_long <- df_long %>%
      group_by(!!!syms(c(feature_id_col, biospecimen_id_col))) %>%
      mutate(n = sum(!is.na(!!sym(measure_col))))
  }
  if (any(df_long$n) > 2) {
    warning("Cannot calculate CV for peptides with 2 or less measurements, removing those peptides")
    df_long <- df_long %>%
      filter(n > 2)
  }

  if ("Step" %in% names(df_long)) {
    if (!is.null(batch_col)) {
      df_long <- df_long %>%
        group_by(!!!syms(c(feature_id_col, batch_col, "Step"))) %>%
        mutate(CV_perBatch = 100 * sd(!!sym(measure_col), na.rm = TRUE) /
          mean(!!sym(measure_col), na.rm = TRUE)) %>%
        ungroup()
    } else {
      warning("batch_col not specified, calculating the total CV only")
    }
    CV_df <- df_long %>%
      group_by(!!!syms(c(feature_id_col, "Step"))) %>%
      mutate(CV_total = 100 * sd(!!sym(measure_col), na.rm = TRUE) /
        mean(!!sym(measure_col), na.rm = TRUE))
    if (!is.null(batch_col)) {
      CV_df <- CV_df %>%
        select(c(!!sym(feature_id_col), CV_total, CV_perBatch)) %>%
        distinct()
    } else {
      CV_df <- CV_df %>%
        select(c(!!sym(feature_id_col), CV_total)) %>%
        distinct()
    }
  } else {
    if (!is.null(batch_col)) {
      df_long <- df_long %>%
        group_by(!!!syms(c(feature_id_col, batch_col))) %>%
        mutate(CV_perBatch = sd(!!sym(measure_col), na.rm = TRUE) /
          mean(!!sym(measure_col), na.rm = TRUE)) %>%
        ungroup()
    } else {
      warning("batch_col not found, calculating the total CV only")
    }
    CV_df <- df_long %>%
      group_by(!!sym(feature_id_col)) %>%
      mutate(CV_total = sd(!!sym(measure_col), na.rm = TRUE) /
        mean(!!sym(measure_col), na.rm = TRUE))
    if (!is.null(batch_col)) {
      CV_df <- CV_df %>%
        select(c(!!sym(feature_id_col), CV_total, CV_perBatch)) %>%
        distinct()
    } else {
      CV_df <- CV_df %>%
        select(c(!!sym(feature_id_col), CV_total)) %>%
        distinct()
    }
  }

  return(CV_df)
}

#' Plot the distribution (boxplots) of per-batch per-step CV of features
#'
#' @inheritParams proBatch
#'
#' @param CV_df data frame with Total CV for each feature & (optionally) per-batch CV
#' @param log_y_scale (logical) whether to display the CV on log-scale
#'
#' @return ggplot object
#' @export
plot_CV_distr.df <- function(CV_df,
                             plot_title = NULL,
                             filename = NULL, theme = "classic", log_y_scale = TRUE) {
  if ("Step" %in% names(CV_df)) {
    gg <- ggplot(CV_df, aes(x = Step, y = CV_total)) +
      geom_boxplot()
  } else {
    gg <- ggplot(CV_df, aes(y = CV_total)) +
      geom_boxplot()
  }
  if (!is.null(plot_title)) {
    gg <- gg + ggtitle(plot_title)
  }
  if (theme == "classic") {
    gg <- gg + theme_classic()
  }
  if (!is.null(filename)) {
    ggsave(gg, filename = filename)
  }

  if (log_y_scale) {
    gg <- gg + scale_y_log10()
  }

  return(gg)
}

#' Plot CV distribution to compare various steps of the analysis
#'
#' @inheritParams proBatch
#' @inheritParams transform_raw_data
#' @param df_long as in \code{df_long} for the rest of the package, but, when it
#' has entries for intensity, represented in \code{measure_col} for several steps,
#' e.g. raw, normalized, batch corrected data, as seen in column \code{Step}, then
#' multi-step CV comparison can be carried out.
#' @param biospecimen_id_col column in \code{sample_annotation}
#' that defines a unique bio ID, which is usually a
#' combination of conditions or groups.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column
#' @param unlog (logical) whether to reverse log transformation of the original data
#'
#' @return \code{ggplot} object with the boxplot of CVs on one or several steps
#' @export
#'
#' @examples
#' data(list = c("example_sample_annotation", "example_proteome"), package = "proBatch")
#' CV_plot <- plot_CV_distr(example_proteome,
#'   sample_annotation = example_sample_annotation,
#'   measure_col = "Intensity", batch_col = "MS_batch",
#'   plot_title = NULL, filename = NULL, theme = "classic"
#' )
plot_CV_distr <- function(df_long, sample_annotation = NULL,
                          feature_id_col = "peptide_group_label",
                          sample_id_col = "FullRunName",
                          measure_col = "Intensity",
                          biospecimen_id_col = "EarTag",
                          batch_col = NULL,
                          unlog = TRUE,
                          log_base = 2,
                          offset = 1,
                          plot_title = NULL,
                          filename = NULL, theme = "classic") {
  CV_df <- calculate_feature_CV(
    df_long = df_long,
    sample_annotation = sample_annotation,
    feature_id_col = feature_id_col,
    sample_id_col = sample_id_col,
    measure_col = measure_col,
    batch_col = batch_col,
    biospecimen_id_col = biospecimen_id_col,
    unlog = unlog,
    log_base = log_base,
    offset = offset
  )
  gg <- plot_CV_distr.df(
    CV_df,
    plot_title = plot_title, filename = filename,
    theme = theme
  )
  return(gg)
}
