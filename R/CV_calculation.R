

#' Calculate CV distribution for each feature
#'
#' @inheritParams proBatch
#'
#' @return data frame with Total CV for each feature & (optionally) per-batch CV
#' @export
#'
#' @examples
#' CV_df = calculate_feature_CV(example_proteome, 
#' sample_annotation = example_sample_annotation, 
#' measure_col = 'Intensity', 
#' batch_col = 'MS_batch')

calculate_feature_CV <- function(df_long, sample_annotation = NULL,
                                 feature_id_col = 'peptide_group_label',
                                 sample_id_col = 'FullRunName',
                                 measure_col = 'Intensity', batch_col = NULL){
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long)
  if (!is.null(batch_col)){
    df_long = df_long %>%
      group_by(!!!syms(c(feature_id_col, batch_col))) %>%
      mutate(n = sum(!is.na(!!sym(measure_col))))
    
  } else {
    df_long = df_long %>%
      group_by(!!sym(feature_id_col)) %>%
      mutate(n = sum(!is.na(!!sym(measure_col))))
  }
  if(any(df_long$n) > 2){
    warning('Cannot calculate CV for peptides with 2 or less measurements, removing those peptides')
    df_long = df_long %>%
      filter(n > 2)
  }
  
  
  if (!is.null(batch_col)){
    df_long = df_long %>%
      group_by(!!!syms(c(feature_id_col, batch_col))) %>%
      mutate(CV_perBatch = sd(!!sym(measure_col), na.rm = T)/mean(!!sym(measure_col), na.rm = T)) %>%
      ungroup()
  } else {
    warning('batch_col not found, calculating the total CV only')
  }
  CV_df = df_long %>%
    group_by(!!sym(feature_id_col)) %>%
    mutate(CV_total = sd(!!sym(measure_col), na.rm = T)/mean(!!sym(measure_col), na.rm = T))
  if(!is.null(batch_col)){
    CV_df = CV_df%>%
      select(c(!!sym(feature_id_col), CV_total, CV_perBatch)) %>%
      distinct()
  } else {
    CV_df = CV_df%>%
      select(c(!!sym(feature_id_col), CV_total)) %>%
      distinct()
  }
  return(CV_df)
}

plot_CV_distr.df <- function(CV_df, 
                          plot_title = NULL, 
                          filename = NULL, theme = 'classic'){
  if ('Step' %in% names(CV_df)){
    gg = ggplot(CV_df,  aes(x = Step, y = CV_total)) +
      geom_boxplot()
  } else {
    gg = ggplot(CV_df,  aes(y = CV_total)) +
      geom_boxplot()
  }
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if(!is.null(filename)){
    ggsave(gg, filename = filename)
  }
  return(gg)
}

#' Plot CV distribution to compare various steps of the analysis
#'
#' @inheritParams proBatch
#' @param df_long as in \code{df_long} for the rest of the package, but, when it 
#' has entries for intensity, represented in \code{measure_col} for several steps, 
#' e.g. raw, normalized, batch corrected data, as seen in column \code{Step}, then
#' multi-step CV comparison can be carried out.
#' @return \code{ggplot} object with the boxplot of CVs on one or several steps
#' @export
#'
#' @examples
#' CV_plot = plot_CV_distr(example_proteome, 
#' sample_annotation = example_sample_annotation, 
#' measure_col = 'Intensity', batch_col = 'MS_batch', 
#' plot_title = NULL, filename = NULL, theme = 'classic')
plot_CV_distr <- function(df_long, sample_annotation = NULL,
                          feature_id_col = 'peptide_group_label',
                          sample_id_col = 'FullRunName',
                          measure_col = 'Intensity', 
                          batch_col = NULL, 
                          plot_title = NULL, 
                          filename = NULL, theme = 'classic'){
  CV_df = calculate_feature_CV(df_long = df_long, 
                               sample_annotation = sample_annotation, 
                               feature_id_col = feature_id_col, 
                               sample_id_col = sample_id_col,
                               measure_col = measure_col, 
                               batch_col = batch_col)
  gg = plot_CV_distr.df(CV_df, plot_title = plot_title, filename = filename, 
                        theme = theme)
  return(gg)
}