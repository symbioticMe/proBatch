

#' Calculate CV distribution for each feature
#'
#' @inheritParams proBatch
#' @param sample_annotation 
#' @param feature_id_col 
#' @param sample_id_col 
#' @param measure_col 
#' @param batch_col 
#'
#' @return data frame with Total CV for each feature & (optionally) per-batch CV
#' @export
#'
#' @examples
#' CV_df = calculate_feature_CV(example_proteome, 
#' sample_annotation = sample_annotation_AgingMice, 
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