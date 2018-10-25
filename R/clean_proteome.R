#' Remove requanted and sparse features
#'
#' Cleans dataset \code{df_long} (\link{proBatch}) by removing requanted features
#' and features not meeting user defined sparsness criterias.
#'
#' @inheritParams proBatch
#' @param missing_frac_batch maximally tolerated fraction of missing values for a
#'   feature in a batch
#' @param missing_frac_total maximally tolerated fraction of globally missing
#'   values for a feature
#'
#' @return \code{df_long} (\link{proBatch}) like data frame filtered as follow:
#'   \itemize{ \item remove requant values \item remove features not meeting
#'   batch or global sparsness thresholds }
#'
#' @export
#'
#' @family dataset cleaning functions
#'
clean_requants <- function(df_long, sample_annotation, peptide_annotation,
                           batch_col = 'MS_batch.final',
                           feature_id_col = 'peptide_group_label',
                           m_score = "m_score",
                           missing_frac_batch = .3, missing_frac_total = .3){
  #for dplyr version 0.7 and higher, this is the way to call the functions
  df_clean = df_long %>%
    merge(peptide_annotation) %>%
    filter(UQ(sym(m_score) != 2) %>%
    merge(sample_annotation) %>%
    group_by_at(vars(one_of(c(c(feature_id_col, batch_col))))) %>%
    mutate(n_samples_in_batch = n()) %>%
    ungroup() %>%
    group_by_at(vars(one_of(batch_col))) %>%
    mutate(n_batch = max(n_samples_in_batch)) %>%
    mutate(requant_fraction_batch = 1 - n_samples_in_batch/n_batch) %>%
    ungroup() %>%
    group_by_at(vars(one_of(c(feature_id_col)))) %>%
    mutate(n_samples_for_peptide = n()) %>%
    mutate(requant_fraction = 1 - n_samples_for_peptide/n_samples) %>%
    filter(requant_fraction < (1 - missing_frac_total)) %>%
    filter(requant_fraction_batch < (1 - missing_frac_batch)) %>%
    ungroup()
  return(df_clean)
}


#' Remove features missing in at least one batch
#'
#' Cleans dataset \code{df_long} (\link{proBatch}) by removing all features that are not present in every batch
#'
#' @inheritParams proBatch
#' @details useful for some downstream functions as ComBat normalization, that
#'   would not work otherwise
#'
#' @return \code{df_long} (\link{proBatch}) like data frame freed of features that were not detected in each batch
#'
#' @export
#'
#' @family dataset cleaning functions
remove_peptides_with_missing_batch <- function(df_long, sample_annotation,
                                               batch_col = 'MS_batch.final',
                                               feature_id_col = 'peptide_group_label'){
  
  n_samples = nrow(sample_annotation)
  features_consistent = df_long %>%
    merge(sample_annotation) %>%
    group_by_at(vars(one_of(c(feature_id_col, batch_col)))) %>%
    summarize(n = n()) %>%
    ungroup () %>%
    complete(!!!syms(c(feature_id_col, batch_col)))%>%
    group_by_at(vars(one_of(c(feature_id_col)))) %>%
    summarize(full_batches = all(!is.na(n))) %>%
    filter(full_batches) %>%
    pull(feature_id_col)

  proteome_clean = df_long %>%
    filter(UQ(sym(feature_id_col)) %in% features_consistent)
  return(proteome_clean)
}


#' Summarize run features
#'
#' Summarizes various peptide properties on a per sample basis. By default will
#' summarize RT, Intensity and m_score. If your feature does not have some of
#' these set them to NULL when calling.
#'
#' @details summarize peptides by sample (ranking) and on the contrary, across
#'   peptide-wise across samples
#'
#' @return a data frame summarizing features in a dataset on a per sample basis.
#'   The following columns are returned: `RT_mean`, `Int_mean`, `numb_requants`,
#'   `median_m_score`, `mean_m_score`, `median_good_m_score` (median of
#'   `m_score` excluding requants)
#'
#' @family dataset cleaning functions
#'
#' @keywords internal
summarize_peptides <- function(df_long, sample_id_col = 'FullRunName',
                               feature_id_col = 'peptide_group_label',
                               RT_col = "RT",
                               measure_col = "Intensity",
                               m_score = "m_score"){
  peptide_summary = df_long %>%
    group_by_at(vars(one_of(sample_id_col)))  %>%
    mutate(rank = rank(measure_col))  %>%
    group_by_at(vars(one_of(feature_id_col))) %>%
    summarise(RT_mean = mean(RT_col ),
              Int_mean = mean(measure_col), rank_mean = mean(rank),
              numb_requants = sum(m_score > 1),
              mean_m_score = mean(m_score),
              median_m_score = median(m_score),
              median_good_m_score = median(m_score[m_score < 1]))
  return(peptide_summary)
}
