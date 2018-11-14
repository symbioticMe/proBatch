cut_into_groups <- function(vect, groups, quantile = T){
  if (quantile){
    cutting_points = unique(quantile(vect, probs = seq(0, 1, length = groups + 1)))
  } else {
    #equally distance RT bins 10 (into number of groups - as argument of the function)
    cut_step = (max(vect)- min(vect))/(groups)
    cutting_points = min(vect)+0:(groups)*cut_step
  }
  res_groups = cut(vect, cutting_points, include.lowest = TRUE)
  return(res_groups)
}

#' sample random peptides for diagnostics
#'
#' @param proteome required columns: \enumerate{
#' \item \code{m_score}
#' \item \code{Intensity}
#' \item \code{peptide_group_label}
#' \item \code{RT}}
#' @param seed
#' @param pep_per_group number of peptides to sample per group
#' @param groups_RT
#' @param groups_intensity
#'
#' @return
#' @export
#'
#' @examples
pick_random_peptides <- function(proteome, seed = 1, pep_per_group = 3,
                                   groups_RT = 10, groups_intensity = 5){
  summarized_proteome = proteome %>%
    filter(m_score < 1) %>%
    mutate(log_intensity = log2(Intensity + 1)) %>%
    group_by(peptide_group_label) %>%
    summarise(RT_average = mean(RT),
              Intensity_average = mean(log_intensity))
  summarized_proteome = summarized_proteome %>%
    mutate(RT_quantile = cut_into_groups(RT_average, groups_RT),
           Int_quantile = cut_into_groups(Intensity_average, groups_intensity)) %>%
    group_by(RT_quantile, Int_quantile)

  set.seed(seed)
  sampled_peptide_df = summarized_proteome %>%
    #sample 3 peptides per each group
    sample_n(pep_per_group) %>%
    arrange(RT_average)
  sampled_peptide_df1 = sampled_peptide_df %>%
    group_by(RT_quantile, Int_quantile) %>% arrange(RT_average) %>%
    mutate(group = rank(peptide_group_label))
  peptides_to_plot = sampled_peptide_df$peptide_group_label
  return(list(peptide_summary_df = summarized_proteome,
         peptides = peptides_to_plot))
}
