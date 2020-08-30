#' Visualise correlation matrix
#'
#' Plot correlation of selected  samples or peptides 
#' @description recommended for heatmap-type visualisation of correlation matrix 
#' with <100 items. With >50 samples and ~10 replicate pairs distribution plots 
#' may be more informative.
#'
#' @inheritParams proBatch
#' @param corr_matrix square correlation matrix
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object
#' @param heatmap_color vector of colors used in heatmap.
#' @param annotation data frame with \code{peptide_annotation} for protein 
#' correlation heatmap or \code{sample_annotation} for sample correlation heatmap
#' @param annotation_id_col \code{feature_id_col} for protein correlation heatmap 
#' or \code{sample_id_col} for sample correlation heatmap
#' @param ... parameters for the \code{\link[pheatmap]{pheatmap}} visualisation,
#'  for details see examples and help to corresponding functions
#'
#' @return \code{pheatmap} object
#' 
#' @export
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}, 
#' \code{\link{plot_sample_corr_distribution}}, 
#' \code{\link{plot_peptide_corr_distribution}}
#' 
#' @examples 
#' peptides <- c("10231_QDVDVWLWQQEGSSK_2", "10768_RLESELDGLR_2")
#' data_matrix_sub = example_proteome_matrix[peptides,]
#' corr_matrix = cor(t(data_matrix_sub), use = 'complete.obs')
#' corr_matrix_plot <- plot_corr_matrix(corr_matrix)
#' 
plot_corr_matrix <- function(corr_matrix,
                             annotation = NULL, 
                             annotation_id_col = 'FullRunName',
                             factors_to_plot = NULL, 
                             cluster_rows = FALSE, cluster_cols = FALSE,
                             heatmap_color = colorRampPalette(
                               rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                             color_list = NULL,
                             filename = NULL, width = 7, height = 7, 
                             units = c('cm','in','mm'),
                             plot_title = NULL, ...) {
  
  #infer the color scheme for annotation (cols & rows)
  if(is.null(color_list) && !is.null(annotation)){
    warning('color_list for annotation (cols & rows) not defined, inferring automatically.
            Numeric/factor columns are guessed, for more controlled color mapping use 
            sample_annotation_to_colors()')
    color_list = sample_annotation_to_colors(sample_annotation = annotation, 
                                             sample_id_col = annotation_id_col, 
                                             factor_columns = factors_to_plot,
                                             numeric_columns = NULL,
                                             guess_factors = TRUE)
  }
  
  if(cluster_rows != cluster_cols){
    warning('different arguments for clustering of rows and columns, this will make
            correlation matrix heatmap assimmetrical!')
  }
  p <- plot_heatmap_generic(corr_matrix, 
                            column_annotation_df = annotation,
                            row_annotation_df = annotation, 
                            fill_the_missing = NULL, 
                            col_ann_id_col = annotation_id_col,
                            row_ann_id_col = annotation_id_col,
                            columns_for_cols = factors_to_plot,
                            columns_for_rows = factors_to_plot,
                            cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                            annotation_color_cols = color_list,
                            annotation_color_rows = color_list,
                            heatmap_color = heatmap_color,
                            filename = filename, width = width, height = height, 
                            units = units, 
                            plot_title = plot_title,
                            ...)
  return(p)
  
}

#' Peptide correlation matrix (heatmap)
#'
#' Plots correlation plot of peptides from a single protein
#'
#' @inheritParams proBatch
#' @param protein_name the name of the protein
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object
#' @param heatmap_color vector of colors used in heatmap.
#' @param ... parameters for the corrplot visualisation
#'
#' @return \code{pheatmap} object
#'
#' @export
#' @examples 
#' protein_corrplot_plot <- plot_protein_corrplot(example_proteome_matrix, 
#' protein_name = 'Haao', peptide_annotation = example_peptide_annotation, 
#' protein_col = 'Gene')
#' 
#' protein_corrplot_plot <- plot_protein_corrplot(example_proteome_matrix, 
#'  protein_name = c('Haao', 'Dhtkd1'), 
#'  peptide_annotation = example_peptide_annotation,
#'  protein_col = 'Gene', factors_to_plot = 'Gene')
#'
plot_protein_corrplot <- function(data_matrix,
                                  protein_name,
                                  peptide_annotation = NULL,
                                  protein_col = 'ProteinName',
                                  feature_id_col = 'peptide_group_label',
                                  factors_to_plot = c('ProteinName'),
                                  cluster_rows = FALSE, cluster_cols = FALSE,
                                  heatmap_color = colorRampPalette(
                                    rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                  color_list = NULL,
                                  filename = NULL,
                                  width = NA, height = NA, 
                                  units = c('cm','in','mm'),
                                  plot_title = sprintf(
                                    'Peptide correlation matrix of %s protein', 
                                    protein_name), ...) {
  
  peptides = peptide_annotation %>%
    filter(!!(sym(feature_id_col)) %in% rownames(data_matrix)) %>%
    filter(!!(sym(protein_col)) %in% protein_name) %>%
    pull(feature_id_col) %>% as.character()
  
  data_matrix_sub = data_matrix[peptides,]
  corr_matrix = cor(t(data_matrix_sub), use = "pairwise.complete.obs")
  
  peptide_annotation = peptide_annotation %>%
    filter(!!(sym(protein_col)) %in% protein_name) %>%
    arrange(!!sym(protein_col))
  
  corr_matrix = corr_matrix[peptide_annotation[[feature_id_col]],
                            peptide_annotation[[feature_id_col]]]
  
  plot_corr_matrix(corr_matrix, 
                   annotation = peptide_annotation, 
                   annotation_id_col = feature_id_col,
                   factors_to_plot = factors_to_plot,
                   cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                   heatmap_color = heatmap_color,
                   color_list = color_list,
                   plot_title = plot_title,
                   filename = filename, width = width, 
                   height = height, units = units, ...)
}

#' Sample correlation matrix (heatmap)
#'
#' Plot correlation of selected samples
#'
#' @inheritParams proBatch
#' @param samples_to_plot string vector of samples in 
#' \code{data_matrix} to be used in the plot
#' @param cluster_rows boolean values determining if rows should be clustered or \code{hclust} object
#' @param cluster_cols boolean values determining if columns should be clustered or \code{hclust} object
#' @param heatmap_color vector of colors used in heatmap.
#' @param ... parameters for the \code{\link[pheatmap]{pheatmap}} visualisation, for details see 
#'   examples and help to corresponding functions
#'
#' @return \code{pheatmap} object
#' 
#' @export
#'
#' @examples
#' specified_samples = example_sample_annotation$FullRunName[
#' which(example_sample_annotation$order %in% 110:115)] 
#' 
#' sample_corr_heatmap <- plot_sample_corr_heatmap(example_proteome_matrix, 
#' samples_to_plot = specified_samples, 
#' factors_to_plot = c('MS_batch','Diet', 'DateTime', 'digestion_batch'),
#'  cluster_rows= FALSE, cluster_cols=FALSE,
#'  annotation_names_col = TRUE, annotation_legend = FALSE, 
#'  show_colnames = FALSE)
#'  
#'  
#'  color_list <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch','EarTag', "Strain", 
#' "Diet", "digestion_batch", "Sex"),
#' numeric_columns = c('DateTime', 'order'))
#'  sample_corr_heatmap_annotated <- plot_sample_corr_heatmap(log_transform_dm(example_proteome_matrix), 
#'  sample_annotation = example_sample_annotation,
#'  factors_to_plot = c('MS_batch','Diet', 'DateTime', 'digestion_batch'),
#'  cluster_rows= FALSE, cluster_cols=FALSE,
#'  annotation_names_col = TRUE, 
#'  show_colnames = FALSE, color_list = color_list)
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}
#' 
plot_sample_corr_heatmap <- function(data_matrix, samples_to_plot = NULL,
                                     sample_annotation = NULL, 
                                     sample_id_col = 'FullRunName',
                                     factors_to_plot = NULL,
                                     cluster_rows = FALSE, cluster_cols = FALSE,
                                     heatmap_color = colorRampPalette(
                                       rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                                     color_list = NULL,
                                     filename = NULL,
                                     width = NA, height = NA, 
                                     units = c('cm','in','mm'),
                                     plot_title = sprintf(
                                       'Correlation matrix of%s samples' ,
                                       ifelse(is.null(samples_to_plot),'',' selected')), ...){
  if(!is.null(samples_to_plot)){
    if (!all(samples_to_plot %in% colnames(data_matrix))){
      missing_samples = setdiff(samples_to_plot, colnames(data_matrix))
      stop(sprintf('The following samples are not in data matrix and can not 
                    be used in sample correlation plotting %s', 
                   paste(missing_samples, collapse = ';\n')))
    }
    corr_matrix = cor(data_matrix[,samples_to_plot], use = 'complete.obs')
  } else {
    corr_matrix = cor(data_matrix, use = 'complete.obs')
  }
  
  if(!is.null(sample_annotation)){
    if (!all(samples_to_plot %in% sample_annotation[[sample_id_col]])){
      warning('some of the samples are not in annotation, this may lead to problems in color annotation')
    }
  }
  
  p <- plot_corr_matrix(corr_matrix, 
                        annotation = sample_annotation, 
                        annotation_id_col = sample_id_col,
                        factors_to_plot = factors_to_plot,
                        cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                        heatmap_color = heatmap_color,
                        color_list = color_list,
                        plot_title = plot_title,
                        filename = filename, width = width, 
                        height = height, units = units, ...)
  return(p)
}

get_sample_corr_df <- function(cor_proteome, sample_annotation,
                               sample_id_col = 'FullRunName',
                               biospecimen_id_col = 'EarTag',
                               batch_col = 'MS_batch'){
  comb_to_keep = data.frame(t(combn(colnames(cor_proteome), 2)))
  names(comb_to_keep) = paste(sample_id_col, seq_len(2), sep = '_')
  
  spec_cols = c(biospecimen_id_col, batch_col)
  
  corr_distribution = melt(cor_proteome,
                           varnames = paste(sample_id_col,seq_len(2), sep = '_'),
                           value.name = 'correlation') %>%
    merge(comb_to_keep) %>%
    merge(sample_annotation %>% select(one_of(c(sample_id_col, spec_cols))),
          by.x = paste(sample_id_col,'1', sep = '_'),
          by.y = sample_id_col, all.x = TRUE) %>%
    data.table::setnames(old = spec_cols, 
                         new = paste(spec_cols, 1, sep = '')) %>%
    merge(sample_annotation %>% select(one_of(c(sample_id_col, spec_cols))),
          by.x = paste(sample_id_col,'2', sep = '_'),
          by.y = sample_id_col, all.x = TRUE) %>%
    data.table::setnames(old = spec_cols, 
                         new = paste(spec_cols, 2, sep = '')) %>%
    mutate(replicate = (!!sym(paste(biospecimen_id_col,'1', sep = '')) == 
                          !!sym(paste(biospecimen_id_col,'2', sep = '')))) %>% 
    mutate(batch_the_same = (!!sym(paste(batch_col,'1', sep = '')) ==
                               !!sym(paste(batch_col,'2', sep = ''))),
           batches = paste(!!sym(paste(batch_col,'1', sep = '')),
                           !!sym(paste(batch_col,'2', sep = '')),sep = ':')) %>%
    mutate(batch_replicate = ifelse(replicate,
                                    ifelse(batch_the_same, 
                                           'same_batch\nsame_biospecimen', 
                                           'same_biospecimen\ndiff_batch'),
                                    ifelse(batch_the_same, 
                                           'same_batch\ndiff_biospecimen',
                                           'diff_batch\ndiff_biospecimen')))
  return(corr_distribution)
}


#' Calculates correlation for all pairs of the samples in data matrix, labels 
#' as replicated/same_batch/unrelated in output columns (see "Value").
#'
#' @inheritParams proBatch
#' @param repeated_samples vector of sample IDs to evaluate, if \code{NULL}, 
#' all samples are taken into account for plotting
#' @param biospecimen_id_col column in \code{sample_annotation} 
#' that defines a unique bio ID, which is usually a 
#' combination of conditions or groups.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column
#'
#' @return dataframe with the following columns, that 
#' are suggested to use for plotting in 
#' \code{\link{plot_sample_corr_distribution}} as \code{plot_param}:
#' \enumerate{
#' \item \code{replicate}
#' \item \code{batch_the_same}
#' \item \code{batch_replicate}
#' \item \code{batches}
#' }
#' other columns are: \enumerate{
#' \item \code{sample_id_1} & \code{sample_id_2}, both 
#' generated from \code{sample_id_col} variable
#' \item \code{correlation} - correlation of two corresponding samples
#' \item \code{batch_1} & \code{batch_2} or analogous, 
#' created the same as \code{sample_id_1}
#' }
#' 
#' @examples 
#' corr_distribution = calculate_sample_corr_distr(data_matrix = example_proteome_matrix, 
#' sample_annotation = example_sample_annotation,
#' batch_col = 'MS_batch',biospecimen_id_col = "EarTag")
#' 
#' @export
#'
calculate_sample_corr_distr <- function(data_matrix, sample_annotation, 
                                        repeated_samples = NULL,
                                        biospecimen_id_col = 'EarTag', 
                                        sample_id_col = 'FullRunName', 
                                        batch_col = 'MS_batch') {
  
  df_long = matrix_to_long(data_matrix, sample_id_col = sample_id_col)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, 
                                     facet_col = NULL, merge = FALSE)
  data_matrix = long_to_matrix(df_long, sample_id_col = sample_id_col)
  
  if (!is.null(repeated_samples)){
    print('calculating correlation of repeated samples only')
    corr_matrix = cor(data_matrix[,repeated_samples], 
                      use = "pairwise.complete.obs")
  } else {
    corr_matrix = cor(data_matrix, use = "pairwise.complete.obs")
  }
  corr_distribution = get_sample_corr_df(cor_proteome = corr_matrix,
                                         sample_annotation = sample_annotation,
                                         sample_id_col = sample_id_col,
                                         biospecimen_id_col = biospecimen_id_col,
                                         batch_col = batch_col)
  return(corr_distribution)
}

#' @name plot_sample_corr_distribution
#' @rdname plot_sample_corr_distribution
#' @title Create violin plot of sample correlation distribution
#'
#' @description Useful to visualize within batch vs within replicate 
#' vs non-related sample correlation
#' 
#' @inheritParams proBatch
#' @param repeated_samples if \code{NULL}, correlation of all samples is plotted
#' @param biospecimen_id_col column in \code{sample_annotation} 
#' that captures the biological sample, 
#' that (possibly) was profiled several times as technical replicates.
#' Tip: if such ID is absent, but can be defined from several columns,
#' create new \code{biospecimen_id} column
#' @param plot_param columns, defined in correlation_df, which is output of
#' \code{calculate_sample_corr_distr}, specifically,  \enumerate{
#' \item \code{replicate}
#' \item \code{batch_the_same}
#' \item \code{batch_replicate}
#' \item \code{batches}
#' }
#' @param corr_distribution data frame with correlation distribution, 
#' as returned by \code{calculate_sample_corr_distr}
#'
#' @return \code{ggplot} type object with violin plot 
#' for each \code{plot_param}
#'
#'

#' 
#' @seealso \code{\link{calculate_sample_corr_distr}}, 
#' \code{\link[ggplot2]{ggplot}}
NULL

#' @rdname plot_sample_corr_distribution
#' 
#' @examples 
#' sample_corr_distribution_plot <- plot_sample_corr_distribution(
#' example_proteome_matrix,
#' example_sample_annotation, batch_col = 'MS_batch', 
#' biospecimen_id_col = "EarTag", 
#' plot_param = 'batch_replicate')
#' 
#' @export
#' 
plot_sample_corr_distribution <- function(data_matrix, sample_annotation,
                                          repeated_samples = NULL,
                                          sample_id_col = 'FullRunName',
                                          batch_col = 'MS_batch',
                                          biospecimen_id_col = 'EarTag',
                                          filename = NULL, width = NA, height = NA, 
                                          units = c('cm','in','mm'),
                                          plot_title = 'Sample correlation distribution',
                                          plot_param = 'batch_replicate',
                                          theme = 'classic'){
    
  if (!is.list(data_matrix)){
    corr_distribution = calculate_sample_corr_distr(data_matrix = data_matrix,
                                                    sample_annotation = sample_annotation,
                                                    sample_id_col = sample_id_col,
                                                    biospecimen_id_col = biospecimen_id_col,
                                                    batch_col = batch_col)
  } else {
    corr_distribution = lapply(seq_len(length(data_matrix)), function(i) {
      dm = data_matrix[[i]]
      corr_distribution = calculate_sample_corr_distr(data_matrix = dm, 
                                                      repeated_samples = repeated_samples, 
                                                      sample_annotation = sample_annotation,
                                                      biospecimen_id_col = biospecimen_id_col,
                                                      sample_id_col =sample_id_col,
                                                      batch_col = batch_col)
      corr_distribution$Step = names(data_matrix)[i]
      return(corr_distribution)
      })
      corr_distribution = do.call(rbind, corr_distribution)
  }
  gg <- plot_sample_corr_distribution.corrDF(corr_distribution = corr_distribution,
                                             filename = filename, 
                                             width = width, height = height, 
                                             units = units,
                                             plot_title = plot_title, 
                                             plot_param = plot_param, 
                                             theme = theme)
  return(gg)
}

#' @rdname plot_sample_corr_distribution
#' 
#' @examples 
#' corr_distribution = calculate_sample_corr_distr(data_matrix = example_proteome_matrix, 
#' sample_annotation = example_sample_annotation,
#' batch_col = 'MS_batch',biospecimen_id_col = "EarTag")
#' sample_corr_distribution_plot <- plot_sample_corr_distribution.corrDF(corr_distribution,
#' plot_param = 'batch_replicate')
#' 
#' \dontrun{
#' sample_corr_distribution_plot <- plot_sample_corr_distribution.corrDF(corr_distribution,
#' plot_param = 'batch_replicate', 
#' filename = 'test_sampleCorr.png', 
#' width = 28, height = 28, units = 'cm')
#' }
#' 
#' @export
#' 
plot_sample_corr_distribution.corrDF <- function(corr_distribution,
                                                 filename = NULL, width = NA, height = NA, 
                                                 units = c('cm','in','mm'),
                                                 plot_title = 'Sample correlation distribution',
                                                 plot_param = 'batch_replicate',
                                                 theme = 'classic', base_size = 20){
  
  gg <- ggplot(corr_distribution, aes_string(x = plot_param, y = 'correlation'))+
    geom_violin(scale = 'width')+
    geom_boxplot(width = .1) +
    theme(axis.title.x=element_blank())
  
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  
  if(('Step' %in% names(corr_distribution)) & 
     length(unique(corr_distribution$Step)) > 1){
    if(length(unique(corr_distribution$Step)) <= 4){
      gg = gg  + facet_grid(.~Step)
    }
    else {
      gg = gg +facet_grid(Step ~ .)
    }
  }
  if (!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic(base_size = base_size)
  } else{
    message("plotting with default ggplot theme, only theme = 'classic' implemented")
  }
  
  if (plot_param =='batches'){
    gg = gg + theme(axis.text.x = element_text(angle = 90))
  }
  if (plot_param == 'batch_replicate'){
    gg = gg + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  }
  
  gg = gg + theme(plot.title = element_text(hjust = .5, face = 'bold'))
  
  save_ggplot(filename, units, width, height, gg)
  return(gg)
}

get_peptide_corr_df <- function(peptide_cor, peptide_annotation, 
                                protein_col = 'ProteinName',
                                feature_id_col = 'peptide_group_label'){
    comb_to_keep = data.frame(t(combn(colnames(peptide_cor), 2)))
    names(comb_to_keep) = paste(feature_id_col, seq_len(2), sep = '_')

    corr_distribution = melt(peptide_cor,
                             varnames = paste(feature_id_col,seq_len(2), sep = '_'),
                             value.name = 'correlation') %>%
      filter(!is.na(correlation)) %>%
      merge(comb_to_keep) %>%
      merge(peptide_annotation %>% 
              select(one_of(c(feature_id_col, protein_col))),
              by.x = paste(feature_id_col,'1', sep = '_'),
              by.y = feature_id_col, all.x = TRUE) %>%
      data.table::setnames(old = protein_col, 
                           new = paste(protein_col, 1, sep = '')) %>%
      merge(peptide_annotation %>% 
              select(one_of(c(feature_id_col, protein_col))),
              by.x = paste(feature_id_col,'2', sep = '_'),
              by.y = feature_id_col, all.x = TRUE) %>%
      data.table::setnames(old = protein_col, 
                           new = paste(protein_col, 2, sep = '')) %>%
      mutate(same_protein = (!!sym(paste(protein_col,'1', sep = '')) ==
                               !!sym(paste(protein_col,'2', sep = '')))) %>%
      mutate(same_protein = ifelse(same_protein, 
                                   'same protein', 'different proteins'))
      
    return(corr_distribution)
}


#' Calculate peptide correlation between and within peptides of one protein
#'
#' @inheritParams proBatch
#'
#' @return dataframe with peptide correlation coefficients 
#' that are suggested to use for plotting in 
#' \code{\link{plot_peptide_corr_distribution}} as \code{plot_param}:
#' 
#' @examples 
#' selected_genes = c('BOVINE_A1ag','BOVINE_FetuinB','Cyfip1')
#' gene_filter = example_peptide_annotation$Gene %in% selected_genes
#' peptides_ann = example_peptide_annotation$peptide_group_label
#' selected_peptides = peptides_ann[gene_filter]
#' matrix_test = example_proteome_matrix[selected_peptides,]
#' pep_annotation_sel = example_peptide_annotation[gene_filter, ]
#' corr_distribution = calculate_peptide_corr_distr(matrix_test, 
#' pep_annotation_sel, protein_col = 'Gene')
#' 
#' @export
#' 
calculate_peptide_corr_distr <- function(data_matrix, peptide_annotation,
                                         protein_col = 'ProteinName',
                                         feature_id_col = 'peptide_group_label'){
    corr_matrix = cor(t(data_matrix), use = "pairwise.complete.obs")
    corr_distribution = get_peptide_corr_df(peptide_cor = corr_matrix,
                                            peptide_annotation = peptide_annotation,
                                            protein_col = protein_col,
                                            feature_id_col = feature_id_col)
    return(corr_distribution)
}

#' @name plot_peptide_corr_distribution
#' @rdname plot_peptide_corr_distribution
#' @title Create violin plot of peptide correlation distribution
#'
#' @description Plot distribution of peptide correlations within one 
#' protein and between proteins
#' 
#' @inheritParams proBatch
#' @param corr_distribution data frame with peptide correlation distribution
#'
#' @return \code{ggplot} object (violin plot of peptide correlation)
#' 
#' @seealso \code{\link{calculate_peptide_corr_distr}}, \code{\link[ggplot2]{ggplot}}
NULL

#' @rdname plot_peptide_corr_distribution
#' 
#' @examples 
#' peptide_corr_distribution <- plot_peptide_corr_distribution(
#' example_proteome_matrix, 
#' example_peptide_annotation, protein_col = 'Gene')
#' 
#' @export
#'
plot_peptide_corr_distribution <- function(data_matrix, peptide_annotation,
                                           protein_col = 'ProteinName',
                                           feature_id_col = 'peptide_group_label',
                                           filename = NULL, width = NA, height = NA, 
                                           units = c('cm','in','mm'),
                                           plot_title = 'Distribution of peptide correlation',
                                           theme = 'classic'){
    
    if (!is.list(data_matrix)){
        corr_distribution = calculate_peptide_corr_distr(data_matrix, 
                                                         peptide_annotation,
                                                         protein_col, 
                                                         feature_id_col)
    } else {
        corr_distribution = lapply(seq_len(length(data_matrix)), function(i) {
            dm = data_matrix[[i]]
            corr_distribution = calculate_peptide_corr_distr(dm, peptide_annotation,
                                                       protein_col, feature_id_col)
            corr_distribution$Step = names(data_matrix)[i]
            return(corr_distribution)
        })
        corr_distribution = do.call(rbind, corr_distribution)%>%
            mutate(Step = factor(Step, levels = names(data_matrix)))
    }
    p = plot_peptide_corr_distribution.corrDF(corr_distribution = corr_distribution,
                                              theme =  theme, 
                                              filename = filename, 
                                              width = width, height = height, 
                                              units = units,
                                              plot_title = plot_title)
    return(p)
}


#' @rdname plot_peptide_corr_distribution
#'
#' @examples 
#' selected_genes = c('BOVINE_A1ag','BOVINE_FetuinB','Cyfip1')
#' gene_filter = example_peptide_annotation$Gene %in% selected_genes
#' peptides_ann = example_peptide_annotation$peptide_group_label
#' selected_peptides = peptides_ann[gene_filter]
#' matrix_test = example_proteome_matrix[selected_peptides,]
#' pep_annotation_sel = example_peptide_annotation[gene_filter, ]
#' corr_distribution = calculate_peptide_corr_distr(matrix_test, 
#' pep_annotation_sel, protein_col = 'Gene')
#' peptide_corr_distribution <- plot_peptide_corr_distribution.corrDF(corr_distribution)
#' 
#' \dontrun{
#' peptide_corr_distribution <- plot_peptide_corr_distribution.corrDF(corr_distribution, 
#' filename = 'test_peptide.png', 
#' width = 28, height = 28, units = 'cm')
#' }
#' 
#' @export
#' 
plot_peptide_corr_distribution.corrDF <- function(corr_distribution, 
                                                  filename = NULL, width = NA, height = NA, 
                                                  units = c('cm','in','mm'),
                                                  plot_title = 'Correlation of peptides', 
                                                  theme = 'classic',
                                                  base_size = 20) {
  median_same_prot = corr_distribution %>%
    filter(same_protein == 'same protein') %>%
    summarize(median = median(correlation, na.rm = TRUE)) %>%
    pull(median)
  gg <- ggplot(corr_distribution, aes_string(x = 'same_protein', 
                                             y = 'correlation'))+
    geom_violin(scale = 'width') +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkgrey') +
    geom_hline(yintercept = median_same_prot, linetype = 'dotted', 
               color = 'tomato1') +
    geom_boxplot(width = .1) +
    xlab(NULL)
  
  if(!is.null(plot_title)){
    gg = gg +
      ggtitle(plot_title)
  }
  
  if(('Step' %in% names(corr_distribution)) & 
     length(unique(corr_distribution$Step)) > 1){
    if(length(unique(corr_distribution$Step)) <= 4){
      gg = gg  + facet_grid(.~Step)
    }
    else {
      gg = gg +facet_grid(Step ~ .)
    }
  }
  
  if (!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic(base_size = base_size)
  } else{
    message("plotting with default ggplot theme, only theme = 'classic' implemented")
  }
  gg = gg + theme(plot.title = element_text(hjust = .5, face = 'bold'))
  
  save_ggplot(filename, units, width, height, gg)
  return(gg)
}
