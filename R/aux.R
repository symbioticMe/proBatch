#' Remove peptides with too many missing peptides
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
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom dplyr one_of
#'
#' @examples
clean_requants <- function(df_long, sample_annotation, batch_column,
                           feature_id_column = 'peptide_group_label',
                           threshold_batch = .7, threshold_global = .5){
  #for dplyr version 0.7 and higher, this is the way to call the functions
  n_samples = nrow(sample_annotation)
  df_clean = df_long %>%
    filter(m_score != 2) %>%
    merge(sample_annotation) %>%
    group_by_at(vars(one_of(c(c(feature_id_column, batch_column))))) %>%
    mutate(n_samples_in_batch = n()) %>%
    ungroup() %>%
    group_by_at(vars(one_of(batch_column))) %>%
    mutate(n_batch = max(n_samples_in_batch)) %>%
    mutate(requant_fraction_batch = 1 - n_samples_in_batch/n_batch) %>%
    ungroup() %>%
    group_by_at(vars(one_of(c(feature_id_column)))) %>%
    mutate(n_samples_for_peptide = n()) %>%
    mutate(requant_fraction = 1 - n_samples_for_peptide/n_samples) %>%
    filter(requant_fraction < threshold_global) %>%
    filter(requant_fraction_batch < threshold_batch) %>%
    ungroup()
  return(df_clean)
}

#' remove peptides that are missing in the whole batch
#' useful for some downstream functions as ComBat normalization, that would "choke"
#'
#' @return data frame free of peptides that were not detected across all batches
#' @export
#' @import dplyr
#' @importFrom tidyr complete
#'
#' @examples
remove_peptides_with_missing_batch <- function(proteome,
                                               batch_column = 'MS_batch.final',
                                               feature_id_column = 'peptide_group_label'){
  peptides_good = proteome %>%
    group_by_at(vars(one_of(c(feature_id_column, batch_column)))) %>%
    summarize(n = n()) %>%
    ungroup () %>%
    complete(!!!rlang::syms(c(feature_id_column, batch_column)))%>%
    group_by_at(vars(one_of(c(feature_id_column)))) %>%
    summarize(full_batches = all(!is.na(n))) %>%
    filter(full_batches) %>%
    pull(feature_id_column)

  proteome_clean = proteome %>%
    filter(rlang::UQ(as.name(feature_id_column)) %in% peptides_good)
  return(proteome_clean)
}

#' Identify stretches of time between runs that are long and split a batches by them
#'
#' @param date_vector
#' @param threshold
#' @param minimal_batch_size
#' @param batch_name
#'
#' @return
#' @export
#'
#' @examples
define_batches_by_MS_pauses <- function(date_vector, threshold,
                                     minimal_batch_size = 5,
                                     batch_name = 'MS_batch'){
  diff = diff(date_vector)
  tipping_points = which(diff > threshold)
  batch_size = diff(tipping_points)
  if (any(batch_size <= minimal_batch_size)){
    warning('some batches are too small, merging with the previous')
    batch_correction = ifelse(batch_size <= minimal_batch_size, batch_size, 0)
    tipping_points = unique(tipping_points+ c(batch_correction, 0))
  }
  tipping_points = c(0, tipping_points, length(date_vector))
  batch_idx = rep(1:(length(tipping_points) -1),
                  times = diff(tipping_points))
  batch_ids = paste(batch_name, batch_idx, sep = '_')
  return(batch_ids)
}


#' summarize peptides by sample (ranking) and on the contrary, across peptide-wise across samples
#'
#' @param proteome
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
summarize_peptides <- function(proteome, sample_id = 'FullRunName',
                               feature_id = 'peptide_group_label'){
  peptide_summary = proteome %>%
    group_by_at(vars(one_of(sample_id)))  %>%
    mutate(rank = rank(Intensity))  %>%
    group_by_at(vars(one_of(feature_id))) %>%
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
#' @import tibble
#' @importFrom magrittr %>%
#' @import dplyr
#' @import reshape2
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
#' @importFrom magrittr %>%
#' @import dplyr
#' @import reshape2
#'
#' @examples
matrix_to_long <- function(data_matrix, sample_annotation, measure_col = 'Intensity', step){
  df_long = data_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = 'peptide_group_label') %>%
    melt(id.var = 'peptide_group_label', value.name = measure_col,
         variable.name = 'FullRunName', factorsAsStrings = F) %>%
    mutate(Step = step) %>%
    merge(sample_annotation)
  return(df_long)
}

#' convert date/time column to POSIX format
#' required to keep number-like behaviour
#'
#' @param sample_annotation
#' @param time_column
#' @param dateTimeFormat
#'
#' @return
#' @export
#'
#' @examples
dates_to_posix <- function(sample_annotation, time_column, new_time_column = NULL,
                          dateTimeFormat = c("%b_%d", "%H:%M:%S")){
  if (length(time_column) == 1){
    if(is.null(new_time_column)) new_time_column = time_column
    time_col = as.character(sample_annotation[[time_column]])
    sample_annotation[[new_time_column]] = as.POSIXct(time_col ,
                                             format=dateTimeFormat)
  }
  else {
    sample_annotation = sample_annotation %>%
      mutate(dateTime = paste(!!!rlang::syms(time_column), sep=" ")) %>%
      mutate(dateTime = as.POSIXct(dateTime,
                                   format = paste(dateTimeFormat, collapse = ' '))) %>%
      rename(!!new_time_column := dateTime)
  }
  return(sample_annotation)
}

#' convert date to order
#'
#' @param sample_annotation
#' @param time_column
#' @param dateTimeFormat
#' @param order_column
#'
#' @return
#' @export
#'
#' @examples
date_to_sample_order <- function(sample_annotation, time_column,
                                 new_time_column = 'DateTime',
                          dateTimeFormat = c("%b_%d", "%H:%M:%S"),
                          order_column = 'order'){
  sample_annotation = dates_to_posix(sample_annotation = sample_annotation,
                                    time_column = time_column,
                                    new_time_column = new_time_column,
                                    dateTimeFormat = dateTimeFormat)
  sample_annotation[[order_column]] = rank(sample_annotation[[new_time_column]])
  return(sample_annotation)
}

#' join list of matrices from different transformation steps into joined data frame
#'
#' @param matrix_list
#' @param Step
#' @param sample_annotation
#' @param measure.col
#'
#' @return
#' @export
#'
#' @examples
join_data_matrices <- function(matrix_list, Step,
                               sample_annotation, measure.col = 'Intensity'){
  long_df_list = lapply(1:length(matrix_list), function(i){
    matrix_to_long(matrix_list[[i]], sample_annotation = sample_annotation,
                   measure_col = measure.col, step = Step[i])
  })
  joined_df = do.call(rbind, long_df_list)

}
