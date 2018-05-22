#' functions for preparing OpenSWATH proteome for batch effect correction and other downstream analyses
#' Used to remove peptides with too many missing values, peptides that are missing in the whole batch (as this causes problems with ComBat) and also for peptide summary
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names) 2) biological and 3) technical covariates (batches etc)
#' @param sample_id_col name of the column in sample_annotation file,
#' where the filenames (colnames of the data matrix are found)
#' @param batch_column column in `sample_annotation` that should be used for batch comparison
#' @param measure_column if `df_long` is among the parameters, it is the column with expression/abundance/intensity,
#' otherwise, it is used internally for consistency
#' @param df_long data frame where each row is a single feature in a single sample,
#' thus it has minimally, `sample_id_col`, `feature_id_column` and `measure_column`, but usually also `m_score` (in OpenSWATH output result file)
#' @param feature_id_column name of the column with feature/gene/peptide/protein ID used with long format matrix (`df_long`). In wide format (`data_matrix`) this would be the row name
#' @name clean_proteome

#' @name clean_proteome
#' @param threshold_batch maximal fraction of missing values for a feature for a batch
#' @param threshold_global maximal fraction of missing values for a feature globally
#'
#' @return `df_long`-like data frame with the requant values removed completely and peptides,
#' missing in batch or globally beyond threshold are removed
#' @export
#' @import dplyr
#' @importFrom magrittr %>%
#' @importFrom dplyr one_of
#'
#' @examples
clean_requants <- function(df_long, sample_annotation,
                           batch_column = 'MS_batch.final',
                           feature_id_column = 'peptide_group_label',
                           threshold_batch = .3, threshold_global = .3){
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
    filter(requant_fraction < (1 - threshold_global)) %>%
    filter(requant_fraction_batch < (1 - threshold_batch)) %>%
    ungroup()
  return(df_clean)
}


#' @details useful for some downstream functions as ComBat normalization, that would not work otherwise
#'
#' @name clean_proteome

#' @return `df_long`-like data frame free of peptides that were not detected in each batch
#' @export
#' @import dplyr
#' @importFrom tidyr complete
#' @examples
remove_peptides_with_missing_batch <- function(df_long,
                                               batch_column = 'MS_batch.final',
                                               feature_id_column = 'peptide_group_label'){
  features_consistent = df_long %>%
    group_by_at(vars(one_of(c(feature_id_column, batch_column)))) %>%
    summarize(n = n()) %>%
    ungroup () %>%
    complete(!!!rlang::syms(c(feature_id_column, batch_column)))%>%
    group_by_at(vars(one_of(c(feature_id_column)))) %>%
    summarize(full_batches = all(!is.na(n))) %>%
    filter(full_batches) %>%
    pull(feature_id_column)

  proteome_clean = df_long %>%
    filter(rlang::UQ(as.name(feature_id_column)) %in% features_consistent)
  return(proteome_clean)
}

#' Identify stretches of time between runs that are long and split a batches by them
#'
#' @param date_vector POSIX or numeric-like vector corresponding to the sample
#' MS profile acquisition timepoint
#' @param threshold time difference that would mean there was an interruption
#' @param minimal_batch_size minimal number of samples in a batch
#' @param batch_name string with a self-explanatory name for the batch (e.g. `MS_batch` for MS-proteomics)
#'
#' @return vector of batches for each sample
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


#' @details summarize peptides by sample (ranking) and on the contrary, across peptide-wise across samples
#'
#' @name clean_proteome
#'
#' @return summarized proteome with columns such as: `RT_mean`, `Int_mean`, `numb_requants`, `median_m_score`, `mean_m_score`, `median_good_m_score` (median of `m_score` other than requants)
#' @export
#' @importFrom magrittr %>%
#'
#' @examples
summarize_peptides <- function(df_long, sample_id_col = 'FullRunName',
                               feature_id_column = 'peptide_group_label'){
  peptide_summary = df_long %>%
    group_by_at(vars(one_of(sample_id_col)))  %>%
    mutate(rank = rank(Intensity))  %>%
    group_by_at(vars(one_of(feature_id_column))) %>%
    summarise(RT_mean = mean(RT),
              Int_mean = mean(Intensity), rank_mean = mean(rank),
              numb_requants = sum(m_score > 1),
              mean_m_score = mean(m_score),
              median_m_score = median(m_score),
              median_good_m_score = median(m_score[m_score < 1]))
  return(peptide_summary)
}
