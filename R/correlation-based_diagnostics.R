#' Visualise correlation matrix
#'
#' Plot correlation of selected  samples or peptides 
#' @description recommended for heatmap-type visualisation of  correlation matrix with
#'  <100 items. With >50 samples and ~10 replicate pairs distribution plots may be more informative.
#'
#' @inheritParams proBatch
#' @param corr_matrix square correlation matrix
#' @param flavor either corrplot from 'corrplot' package or 
#' heatmap, as in 'pheatmap'
#' @param filename path where the results are saved. 
#' If null the object is returned to the active window;
#' otherwise, the object is save into the file. Currently only 
#' pdf and png format is supported
#' @param width option  determining the output image width
#' @param height option  determining the output image width
#' @param unit units: 'cm', 'in' or 'mm'
#' @param plot_title Title of the correlation plot (e.g., processing step + 
#' representation level (fragments, transitions, proteins) or Protein Name)
#' @param ... parameters for the \code{\link[corrplot]{corrplot.mixed}} or
#' \code{\link[pheatmap]{pheatmap}} visualisation, for details see examples and
#'   help to corresponding functions
#'
#' @return \code{corrplot} or \code{pheatmap} object depending on \code{flavor}
#' 
#' @export
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}, \code{\link[corrplot]{corrplot.mixed}},
#' \code{\link{plot_sample_corr_distribution}}, \code{\link{plot_peptide_corr_distribution}}
#' 
#' @examples 
#' peptides <- c("10231_QDVDVWLWQQEGSSK_2", "10768_RLESELDGLR_2")
#' data_matrix_sub = example_proteome_matrix[peptides,]
#' corr_matrix = cor(t(data_matrix_sub), use = 'complete.obs')
#' corr_matrix_plot <- plot_corr_matrix(corr_matrix,  flavor = "corrplot")
#' 
plot_corr_matrix <- function(corr_matrix, flavor = c('pheatmap','corrplot'), 
                             filename = NULL, width = NA, height = NA, 
                             unit = c('cm','in','mm'),
                             plot_title = NULL, ...) {
  
  flavor <- match.arg(flavor)    
  
  if (!(flavor %in% c('pheatmap','corrplot'))){
    stop('only pheatmap or corrplot can be produced to illustrate sample correlation,
           choose one of these two options')
  }
  if (is.null(filename)){
    switch(flavor,
           corrplot = corrplot.mixed(corr_matrix, title = plot_title, ...),
           pheatmap = pheatmap(corr_matrix, main = plot_title, ...))
  } else {
    if(length(unit) > 1) unit = unit[1]
    if (unit == 'mm'){
      unit = 'cm'
      width = width/10
      height = height/10
    }
    if(unit == 'cm'){
      width  = width/2.54
      height = height/2.54
    }
    if (flavor == 'corrplot'){
      pdf(file = filename, width = width, height = height, title = plot_title)
      corrplot.mixed(corr_matrix, title = plot_title, ...)
      dev.off()
    } else{
      pheatmap(corr_matrix, filename = paste(filename, '.pdf', sep = ''),
               width = width, height = height, main = plot_title, ...)
    }
    
  }
}

#' Peptide correlation matrix (heatmap)
#'
#' Plots correlation plot of peptides from a single protein
#'
#' @inheritParams proBatch
#' @param protein_name the name of the protein
#' @param flavor either corrplot from 'corrplot' 
#' package or heatmap, as in 'pheatmap'
#' @param filename path where the results are saved. 
#' If null the object is returned to the active window;
#' otherwise, the object is save into the file. Currently 
#' only pdf and png format is supported
#' @param width option  determining the output image width
#' @param height option  determining the output image width
#' @param unit units: 'cm', 'in' or 'mm'
#' @param plot_title The title of the plot, e.g. protein name and the processing step
#' @param ... parameters for the corrplot visualisation
#'
#' @return \code{corrplot} or \code{pheatmap} object depending on \code{flavor}
#'
#' @export
#' @examples 
#' protein_corrplot_plot <- plot_protein_corrplot(example_proteome_matrix, protein_name = 'Haao',
#'                peptide_annotation = example_peptide_annotation, 
#'                protein_col = 'Gene', flavor = "pheatmap")
#'
plot_protein_corrplot <- function(data_matrix,
                                  protein_name,
                                  peptide_annotation,
                                  protein_col = 'ProteinName',
                                  feature_id_col = 'peptide_group_label',
                                  flavor = c('pheatmap','corrplot'),
                                  filename = NULL,
                                  width = NA, height = NA, unit = c('cm','in','mm'),
                                  plot_title = sprintf(
                                    'Peptide correlation matrix of %s protein', 
                                    protein_name), ...) {
  
  flavor <- match.arg(flavor)    
  peptides = peptide_annotation %>%
    filter(!!(sym(feature_id_col)) %in% rownames(data_matrix)) %>%
    filter(!!(sym(protein_col)) == protein_name) %>%
    pull(feature_id_col) %>% as.character()
  data_matrix_sub = data_matrix[peptides,]
  corr_matrix = cor(t(data_matrix_sub), use = 'complete.obs')
  plot_corr_matrix(corr_matrix, plot_title = plot_title, flavor = flavor,
                   filename = filename, width = width, 
                   height = height, unit = unit, ...)
}

#' Sample correlation matrix (heatmap)
#'
#' Plot correlation of selected samples
#'
#' @inheritParams proBatch
#' @param samples_to_plot string vector of samples in 
#' \code{data_matrix} to be used in the plot
#' @param filename path where the results are saved. 
#' If null the object is returned to the active window;
#' otherwise, the object is save into the file. 
#' Currently only pdf and png format is supported
#' @param width option  determining the output image width
#' @param height option  determining the output image width
#' @param flavor either corrplot from 'corrplot' package or 
#' heatmap, as in 'pheatmap'
#' @param unit units: 'cm', 'in' or 'mm'
#' @param plot_title Title of the plot (usually, processing step + 
#' representation level (fragments, transitions, proteins))
#' @param ... parameters for the \code{\link[corrplot]{corrplot.mixed}} or
#' \code{\link[pheatmap]{pheatmap}} visualisation, for details see 
#'   examples and help to corresponding functions
#'
#' @return \code{corrplot} or \code{pheatmap} object depending on \code{flavor}
#' 
#' @export
#'
#' @examples
#' specified_samples = example_sample_annotation$FullRunName[
#' which(example_sample_annotation$order %in% 110:115)] 
#' 
#' sample_corr_heatmap <- plot_sample_corr_heatmap(example_proteome_matrix, 
#' samples_to_plot = specified_samples, 
#'  flavor = 'pheatmap', 
#'  cluster_rows= FALSE, cluster_cols=FALSE,
#'  annotation_names_col = TRUE, annotation_legend = FALSE, 
#'  show_colnames = FALSE)
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}, 
#' \code{\link[corrplot]{corrplot.mixed}}
#' 
plot_sample_corr_heatmap <- function(data_matrix, samples_to_plot = NULL,
                                     flavor = c('pheatmap','corrplot'), filename = NULL,
                                     width = NA, height = NA, 
                                     unit = c('cm','in','mm'),
                                     plot_title = sprintf(
                                       'Correlation matrix of sample %s',
                                       samples_to_plot), ...){
  flavor <- match.arg(flavor)    
  if(!is.null(samples_to_plot)){
    corr_matrix = cor(data_matrix[,samples_to_plot], use = 'complete.obs')
  } else {
    corr_matrix = cor(data_matrix, use = 'complete.obs')
  }
  plot_corr_matrix(corr_matrix, plot_title = plot_title, flavor = flavor,
                   filename = filename, width = width, 
                   height = height, unit = unit, ...)
}


#' Transform square correlation matrix into long data.frame of correlations
#' of the replicated samples, samples from same batch and unrelated samples
#' @description useful to infer summary statistics (medians etc.) of correlation 
#' within/between replicates and batches.
#'
#' @inheritParams proBatch
#' @param cor_proteome sample correlation matrix (square)
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
#' @export
#' 
get_sample_corr_df <- function(cor_proteome, sample_annotation,
                                    sample_id_col = 'FullRunName',
                                    biospecimen_id_col = 'EarTag',
                                    batch_col = 'batch'){
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
    data.table::setnames(old = spec_cols, new = paste(spec_cols, 1, sep = '')) %>%
    merge(sample_annotation %>% select(one_of(c(sample_id_col, spec_cols))),
          by.x = paste(sample_id_col,'2', sep = '_'),
          by.y = sample_id_col, all.x = TRUE) %>%
    data.table::setnames(old = spec_cols, new = paste(spec_cols, 2, sep = '')) %>%
    mutate(replicate = (!!sym(paste(biospecimen_id_col,'1', sep = '')) == 
                          !!sym(paste(biospecimen_id_col,'2', sep = '')))) %>% 
    mutate(batch_the_same = (!!sym(paste(batch_col,'1', sep = '')) ==
                               !!sym(paste(batch_col,'2', sep = ''))),
           batches = paste(!!sym(paste(batch_col,'1', sep = '')),
                           !!sym(paste(batch_col,'2', sep = '')), sep = ':')) %>%
    mutate(batch_replicate = ifelse(replicate,
                                    ifelse(batch_the_same, 'same_batch\nsame_biospecimen', 
                                           'same_biospecimen\ndiff_batch'),
                                    ifelse(batch_the_same, 'same_batch\ndiff_biospecimen',
                                           'diff_batch\ndiff_biospecimen')))
  return(corr_distribution)
}


#' Calculates correlation of data matrix and calculates correlation distribution for all pairs 
#' of the replicated samples 
#'
#' @inheritParams proBatch
#' @param repeated_samples if \code{NULL}, only repeated sample correlation is plotted
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
#' @export
#'
calculate_sample_corr_distribution <- function(data_matrix, repeated_samples, sample_annotation,
                              biospecimen_id_col, sample_id_col, batch_col) {
  if (!is.null(repeated_samples)){
    print('calculating correlation of repeated samples only')
    corr_matrix = cor(data_matrix[,repeated_samples], use = "pairwise.complete.obs")
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

#' Create violin plot of correlation distribution
#'
#' Useful to visualize within batch vs within replicate 
#' vs non-related sample correlation
#' 
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data matrix with 1) \code{sample_id_col} (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param repeated_samples if \code{NULL}, only repeated sample correlation is plotted
#' @param biospecimen_id_col column in \code{sample_annotation} 
#' that captures the biological sample, 
#' that (possibly) was profiled several times as technical replicates.
#' Tip: if such ID is absent, but can be defined from several columns,
#' create new \code{biospecimen_id} column
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix) are found
#' @param batch_col column in \code{sample_annotation} that should be used for
#'   batch comparison
#' @param plot_param columns, defined in correlation_df, which is output of
#' \code{get_sample_corr_distrib}, specifically, #' \enumerate{
#' \item \code{replicate}
#' \item \code{batch_the_same}
#' \item \code{batch_replicate}
#' \item \code{batches}
#' };
#'
#' @return \code{ggplot} type object with violin plot 
#' for each \code{plot_param}
#'
#' @export
#'
#' @examples 
#' sample_corr_distribution_plot <- plot_sample_corr_distribution(
#' example_proteome_matrix,
#' example_sample_annotation, batch_col = 'MS_batch', 
#' biospecimen_id_col = "EarTag", 
#' plot_param = 'batch_replicate')
#' 
#' @seealso \code{\link{get_sample_corr_distrib}}, \code{\link[ggplot2]{ggplot}}
plot_sample_corr_distribution <- function(data_matrix, sample_annotation,
                                          repeated_samples = NULL,
                                          sample_id_col = 'FullRunName',
                                          batch_col = 'MS_batch',
                                          biospecimen_id_col = 'EarTag',
                                          plot_title = 'Sample correlation distribution',
                                          plot_param = 'batch_replicate',
                                          theme = 'classic'){
    
    if(!setequal(unique(sample_annotation[[sample_id_col]]), 
                unique(colnames(data_matrix)))){
        warning('Sample IDs in sample annotation not 
                consistent with samples in input data.')}
    
    if (!is.list(data_matrix)){
        corr_distribution = calculate_sample_corr_distribution(data_matrix = data_matrix, 
                                              repeated_samples = repeated_samples, 
                                              sample_annotation = sample_annotation,
                                              biospecimen_id_col = biospecimen_id_col, 
                                              sample_id_col = sample_id_col, 
                                              batch_col = batch_col)
    } else {
        corr_distribution = lapply(1:length(data_matrix), function(i) {
            dm = data_matrix[[i]]
            corr_distribution = calculate_sample_corr_distribution(data_matrix = dm, 
                                                  repeated_samples = repeated_samples, 
                                                  sample_annotation = sample_annotation,
                                                  biospecimen_id_col = biospecimen_id_col, 
                                                  sample_id_col =sample_id_col, batch_col = batch_col)
            corr_distribution$Step = names(data_matrix)[i]
            return(corr_distribution)
        })
        corr_distribution = do.call(rbind, corr_distribution)
    }
    p <- plot_sample_corr_distribution.corrDF(corr_distribution = corr_distribution,
                                              plot_title = plot_title, plot_param = plot_param, 
                                              theme = theme)
    return(p)
}

#' Plot correlation distribution, starting from the correlation matrix
#'
#' @inheritParams proBatch
#' @param corr_distribution data frame with correlation distribution
#' @param plot_param one of the columns, as returned by \code{get_sample_corr_df}
#'
#' @return \code{ggplot} type object with violin plot
#' @export
#'
#' @examples
plot_sample_corr_distribution.corrDF <- function(corr_distribution,
                                                 plot_title = 'Sample correlation distribution',
                                                 plot_param = 'batch_replicate',
                                                 theme = 'classic'){
  
  gg <- ggplot(corr_distribution, aes_string(x = plot_param, y = 'correlation'))+
    geom_violin(scale = 'width')+
    geom_boxplot(width = .1) 
  
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  
  if(('Step' %in% names(corr_distribution)) & length(unique(corr_distribution$Step)) > 1){
    if(length(unique(corr_distribution$Step)) <= 4){
      gg = gg  + facet_grid(.~Step)
    }
    else {
      gg = gg +facet_grid(Step ~ .)
    }
  }
  if (theme == 'classic'){
    gg = gg + theme_classic()
  }
  
  if (plot_param =='batches'){
    gg = gg + theme(axis.text.x = element_text(angle = 90))
  }
  if (plot_param == 'batch_replicate'){
    gg = gg + theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  }
  
  gg = gg + theme(plot.title = element_text(hjust = .5, face = 'bold'))
  return(gg)
}



#' Transform square correlation matrix into long data frame of correlations
#'
#' @inheritParams proBatch
#' @param peptide_cor peptide correlation matrix (square)
#'
#' @return dataframe with correlation coefficients for each peptide pair and label
#' \code{same_protein} or \code{different_proteins} in column \code{same_protein}. 
#' Intended for use in \code{\link{plot_peptide_corr_distribution}}.
#' 
#' @export
#' 
get_peptide_corr_df <- function(peptide_cor, peptide_annotation, protein_col = 'ProteinName',
                                feature_id_col = 'peptide_group_label'){
    comb_to_keep = data.frame(t(combn(colnames(peptide_cor), 2)))
    names(comb_to_keep) = paste(feature_id_col, 1:2, sep = '_')

    corr_distribution = melt(peptide_cor,
                             varnames = paste(feature_id_col,1:2, sep = '_'),
                             value.name = 'correlation') %>%
      merge(comb_to_keep) %>%
      merge(peptide_annotation %>% select(one_of(c(feature_id_col, protein_col))),
              by.x = paste(feature_id_col,'1', sep = '_'),
              by.y = feature_id_col, all.x = TRUE) %>%
      data.table::setnames(old = protein_col, new = paste(protein_col, 1, sep = '')) %>%
      merge(peptide_annotation %>% select(one_of(c(feature_id_col, protein_col))),
              by.x = paste(feature_id_col,'2', sep = '_'),
              by.y = feature_id_col, all.x = TRUE) %>%
      data.table::setnames(old = protein_col, new = paste(protein_col, 2, sep = '')) %>%
      mutate(same_protein = (!!sym(paste(protein_col,'1', sep = '')) ==
                               !!sym(paste(protein_col,'2', sep = '')))) %>%
      mutate(same_protein = ifelse(same_protein, 'same protein', 'different proteins'))
      
    return(corr_distribution)
}


#' Transform square correlation matrix into long data frame of correlations
#'
#' @inheritParams proBatch
#'
#' @return dataframe with peptide correlation coefficients 
#' that are suggested to use for plotting in 
#' \code{\link{plot_peptide_corr_distribution}} as \code{plot_param}:
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

#' Plot distribution of peptide correlations within one 
#' protein and between proteins
#'
#' @inheritParams proBatch
#'
#' @return \code{ggplot} type object with violin plot
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
                                           plot_title = 'Distribution of peptide correlation',
                                           theme = 'classic'){
    
    if (!is.list(data_matrix)){
        corr_distribution = calculate_peptide_corr_distr(data_matrix, peptide_annotation,
                                                   protein_col, feature_id_col)
    } else {
        corr_distribution = lapply(1:length(data_matrix), function(i) {
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
                                              theme =  theme, plot_title = plot_title)
    return(p)
}

#' Plot distribution of peptide correlation within and between proteins, 
#' starting from peptide correlation data frame
#'
#' @inheritParams proBatch
#' @param corr_distribution 
#'
#' @return \code{ggplot} type object with violin plot
#' 
#' @export
#'
#' @examples
plot_peptide_corr_distribution.corrDF <- function(corr_distribution, 
                                                  plot_title = 'Correlation of peptides', 
                                                  theme = 'classic') {
  median_same_prot = corr_distribution %>%
    filter(same_protein == 'same protein') %>%
    summarize(median = median(correlation, na.rm = T)) %>%
    pull(median)
  gg <- ggplot(corr_distribution, aes_string(x = 'same_protein', y = 'correlation'))+
    geom_violin(scale = 'width') +
    geom_hline(yintercept = 0, linetype = 'dashed', color = 'darkgrey') +
    geom_hline(yintercept = median_same_prot, linetype = 'dotted', color = 'tomato1') +
    geom_boxplot(width = .1) +
    xlab(NULL)
  
  if(!is.null(plot_title)){
    gg = gg +
      ggtitle(plot_title)
  }
  
  if(('Step' %in% names(corr_distribution)) & length(unique(corr_distribution$Step)) > 1){
    if(length(unique(corr_distribution$Step)) <= 4){
      gg = gg  + facet_grid(.~Step)
    }
    else {
      gg = gg +facet_grid(Step ~ .)
    }
  }
  
  if (theme == 'classic'){
    gg = gg + theme_classic()
  }
  gg = gg + theme(plot.title = element_text(hjust = .5, face = 'bold'))
  return(gg)
}
