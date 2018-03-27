library(WGCNA)

#' Title
#'
#' @param df_long
#' @param sample_annotation
#' @param batch_column
#' @param feature_id_column
#' @param threshold_batch
#' @param threshold_global
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
clean_requants <- function(df_long, sample_annotation, batch_column,
                           feature_id_column = 'peptide_group_label',
                           threshold_batch = .7, threshold_global = .5){
  #for dplyr version 0.7 and higher, this is the way to call the functions
  feature_id_column_var <- quo(feature_id_column)

  peptide_summary = df_long %>%
    merge(sample_annotation) %>%
    dplyr::select(one_of(c(feature_id_column, batch_column, 'm_score'))) %>%
    group_by(!!(c(feature_id_column, batch_column))) %>%
    mutate(n_requants_batch = sum(m_score == 2), n_batch = n()) %>%
    mutate(requant_fraction_batch = n_requants_batch/n_batch) %>%
    group_by_at(vars(one_of(c(feature_id_column)))) %>%
    mutate(n_requants = sum(m_score == 2), n = n()) %>%
    mutate(requant_fraction = n_requant/n) %>%
    filter(requant_fraction < threshold_global) %>%
    filter(requant_fraction_batch < threshold_batch)
  return(df_long %>% filter(!!feature_id_column %in% peptide_summary[[feature_id_column]]))
}

sample_annotation_to_colors <- function(sample_annotation, columns_to_include = NULL,
                                        columns_to_exclude = NULL){

}

#' Title
#'
#' @param proteome
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
summarize_peptides <- function(proteome){
  peptide_summary = proteome %>%
    group_by(FullRunName) %>%
    mutate(rank = rank(Intensity))  %>%
    group_by(peptide_group_label, ProteinName) %>%
    summarise(RT_mean = mean(RT),
              Int_mean = mean(Intensity), rank_mean = mean(rank),
              numb_requants = sum(m_score > 1),
              mean_m_score = mean(m_score),
              median_m_score = median(m_score),
              median_good_m_score = median(m_score[m_score < 1]))
}

#' Convert from long data frame to data matrix (features in rows, samples in columns)
#'
#' @param proteome_long
#' @param feature_id_column
#' @param measure_column
#' @param sample_id_column
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
convert_to_matrix <- function(data_df_long,
                              feature_id_column = 'peptide_group_label',
                              measure_column = 'Intensity',
                              sample_id_column = 'FullRunName'){
  casting_formula =  as.formula(paste(feature_id_column, sample_id_column,
                                      sep =  " ~ "))
  proteome_wide = dcast(data_df_long, formula=casting_formula,
                        value.var=measure_column) %>%
    column_to_rownames(feature_id_column) %>%
    as.matrix()
  return(proteome_wide)
}

#' Convert the features x samples data matrix to a long format (e.g. for plotting)
#'
#' @param data_matrix
#' @param step
#' @param sample_annotation
#'
#' @return
#' @export
#' @import tidyverse, reshape2
#'
#' @examples
matrix_to_long <- function(data_matrix, sample_annotation, measure_col = 'Intensity', step){
  df_long = data_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = 'peptide_group_label') %>%
    melt(id.var = 'peptide_group_label', value.name = measure_col,
         variable.name = 'FullRunName', factorsAsStrings = F) %>%
    mutate(before_after = step) %>%
    merge(sample_annotation)
  return(df_long)
}

#' Add order and dateTime column to sample annotation data frame
#'
#' @param sample_annotation
#' @param time_column
#' @param dateTimeFormat
#' @param order_column
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
date_to_order <- function(sample_annotation, time_column,
                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                          order_column = 'order'){
  if (length(time_column) == 1){
    sample_annotation$run_order = as.POSIXct(sample_annotation[[time_column]],
                                             format=dateTimeFormat)
  }
  else {
    sample_annotation = sample_annotation %>%
      mutate(dateTime = as.POSIXct(paste(), paste(dateTimeFormat)))
  }

  #sample_annotation$run_order = as.POSIXct(paste(sample_annotation$RunDate,
  #                                               sample_annotation$RunTime),
  #                                         format=dateTimeFormat)
  sample_annotation[[order_column]] = rank(sample_annotation$run_order)
  return(sample_annotation)
}
