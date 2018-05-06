#' plot correlation plot of a single protein
#'
#' @param data_matrix
#' @param peptide_name
#' @param title
#'
#' @return
#' @export
#' @import corrplot
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples plot_corr_plot(q_norm_proteome, protein_name = 'Haao',
#' peptide_annotation = peptide_annotation, prot.column = 'Gene',
#' title = 'Haao protein peptides after quantile norm',
#' number.cex=0.75, tl.cex = .75
#' mar=c(0,0,1,0))
plot_corr_plot_protein <- function(data_matrix, protein_name, peptide_annotation,
                           prot.column = 'ProteinName',
                           peptide_col_name = 'peptide_group_label',
                           title = NULL, ...){
  peptides = peptide_annotation %>%
    filter(rlang::UQ(as.name(peptide_col_name)) %in% rownames(data_matrix)) %>%
    filter(rlang::UQ(as.name(prot.column)) == protein_name) %>%
    pull(peptide_col_name)
    #peptide_annotation[[peptide_col_name]][peptide_annotation[[prot.column]] == protein_name]
  data_matrix_sub = data_matrix[peptides, ]
  corr_matrix = within_prot_correlation(data_matrix_sub)
  corrplot.mixed(corr_matrix, lower = "ellipse", upper = "number",
                 tl.col = "black", diag = 'l', tl.pos = "lt", ...)
}

#' Plot correlation of selected samples
#'
#' @param data_matrix features x samples matrix, with sample IDs as colnames
#' @param samples_to_plot string vector of samples from data_matrix
#' @param ... parameters for the corrplot visualisation
#'
#' @return
#' @export
#' @importFrom corrplot corrplot.mixed
#' @importFrom pheatmap pheatmap
#'
#' @examples
plot_corr_between_samples <- function(data_matrix, samples_to_plot,
                                      flavor = 'corrplot', ...){
  corr_matrix = cor(data_matrix[,samples_to_plot], use = 'complete.obs')
  switch(flavor,
         corrplot = corrplot.mixed(corr_matrix, lower = "ellipse", upper = "number",
                                   tl.col = "black", diag = 'l', tl.pos = "lt",
                                   number.cex=0.75, tl.cex = .75, ...),
         pheatmap = pheatmap(corr_matrix, ...))


}


#' calculate correlation distribution for all pairs of the replicated samples
#'
#' @param cor_proteome
#' @param sample_annotation
#' @param sample_id_col
#' @param biospecimen_id_col
#' @param batch_col
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#'
#'
#' @examples
get_sample_corr_distrib <- function(cor_proteome, sample_annotation,
                                   sample_id_col = 'FullRunName',
                                   biospecimen_id_col = 'EarTag',
                                   batch_col = 'MS_batch.final'){

  comb_to_keep = data.frame(t(combn(colnames(cor_proteome), 2)))
  names(comb_to_keep) = paste(sample_id_col, 1:2, sep = '_')

  spec_cols = c(biospecimen_id_col, batch_col)

  corr_distribution = melt(cor_proteome,
                           varnames = paste(sample_id_col,1:2, sep = '_'),
                           value.name = 'correlation') %>%
    merge(comb_to_keep) %>%
    merge(sample_annotation %>% select(one_of(c(sample_id_col, spec_cols))),
          by.x = paste(sample_id_col,'1', sep = '_'),
          by.y = sample_id_col, all.x = T) %>%
    data.table::setnames(old = spec_cols, new = paste(spec_cols, 1, sep = '')) %>%
    merge(sample_annotation %>% select(one_of(c(sample_id_col, spec_cols))),
          by.x = paste(sample_id_col,'2', sep = '_'),
          by.y = sample_id_col, all.x = T) %>%
    data.table::setnames(old = spec_cols, new = paste(spec_cols, 2, sep = '')) %>%
    mutate(replicate = (!!rlang::sym(paste(biospecimen_id_col,'1', sep = '')) ==
                          !!rlang::sym(paste(biospecimen_id_col,'2', sep = '')))) %>%
    mutate(batch_the_same = (!!rlang::sym(paste(batch_col,'1', sep = '')) ==
                               !!rlang::sym(paste(batch_col,'2', sep = ''))),
           batches = paste(!!rlang::sym(paste(batch_col,'1', sep = '')),
                           !!rlang::sym(paste(batch_col,'2', sep = '')), sep = ':')) %>%
    mutate(batch_replicate = ifelse(replicate,
                                    ifelse(batch_the_same, 'same_batch_same_replicate', 'same_replicate_diff_batch'),
                                    ifelse(batch_the_same, 'same_batch_diff_biospecimen','diff_batch_diff_biospecimen')))
  return(corr_distribution)
}

#' Create violin plot of correlation distribution
#' Typically to visualize within batch vs within replicate vs non-related sample correlation
#'
#' @param repeated_samples
#' @param covariate
#' @param title
#' @param sample_id_col
#' @param batch_col
#' @param plot_param
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
plot_sample_corr_distribution <- function(data_matrix, sample_annotation,
                                   repeated_samples = NULL,
                                   sample_id_col = 'FullRunName',
                                   batch_col = 'MS_batch.final',
                                   covariate = 'EarTag',
                                   title = 'Correlation_distribution',
                                   plot_param = 'batch_the_same'){
  corr_distribution <- function(data_matrix, repeated_samples, sample_annotation,
                                covariate) {
    corr_matrix = cor(data_matrix[,repeated_samples], use = 'complete.obs')
    corr_distribution = get_sample_corr_distrib(cor_proteome = corr_matrix,
                                                sample_annotation = sample_annotation,
                                                sample_id_col = sample_id_col,
                                                biospecimen_id_col = covariate,
                                                batch_col = batch_col)
    return(corr_distribution)
  }
  if (!is.list(data_matrix)){
    corr_distribution = corr_distribution(data_matrix, repeated_samples, sample_annotation, covariate)
  } else {
    corr_distribution = lapply(1:length(data_matrix), function(i) {
      dm = data_matrix[[i]]
      corr_distribution = corr_distribution(dm, repeated_samples, sample_annotation, covariate)
      corr_distribution$Step = names(data_matrix)[i]
      return(corr_distribution)
    })
    corr_distribution = do.call(rbind, corr_distribution)
  }
  p <- ggplot(corr_distribution, aes_string(x = plot_param, y = 'correlation'))+
    geom_violin()+
    geom_boxplot(width = .1)+
    ggtitle(title)+
    theme_classic()
  if (plot_param =='batches'){
    p = p + theme(axis.text.x = element_text(angle = 90))
  }
  if(is.list(data_matrix)){
    if(length(data_matrix) < 4){
      p = p + theme_bw() + facet_grid(.~Step)
    }
    else {
      p = p + theme_bw()+facet_grid(Step ~ .)
    }
  }
  p = p + theme(plot.title = element_text(hjust = .5, face = 'bold'))
  return(p)
}



within_prot_correlation <- function(data_matrix_sub){
  cor_matrix = cor(t(data_matrix_sub), use = 'complete.obs')
  return(cor_matrix)
}

get_prot_corr_df <- function(data_martix){
  cor_matrix = within_prot_correlation(data_martix)
  cor_vector = cor_matrix[upper.tri(cor_matrix)]
  corr_df = data.frame(correlation = cor_vector)
  return(corr_df)
}

#' Plot distribution of correlations
#'
#' @param data_matrix_sub
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
distribution_of_cor <- function(data_matrix_sub, facet_var = NULL, theme = 'classic'){
  corr_df = get_prot_corr_df(data_matrix_sub)
  gg = ggplot(corr_df,
              aes(x = correlation))+
    geom_histogram()
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if(!is.null(facet_var)){
    if (length(facet_var) == 1){
      gg = gg + facet_wrap(as.formula(paste("~", facet_var)))
    } else {
      if(length(facet_var) == 2){
        gg = gg + facet_grid(as.formula(paste(facet_var, collapse = ' ~ ')))
      }
    }

  }
  return(gg)
}
