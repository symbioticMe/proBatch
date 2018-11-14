#' Visualise correlation matrix
#'
#' Plot correlation of selected  samples or peptides
#'
#' @inheritParams proBatch
#' @param corr_matrix square correlation matrix
#' @param flavor either corrplot from 'corrplot' package or heatmap, as in 'pheatmap'
#' @param filename path where the results are saved. If null the object is returned to the active window;
#' otherwise, the object is save into the file. Currently only pdf and png format is supported
#' @param width option  determining the output image width
#' @param height option  determining the output image width
#' @param unit units: 'cm', 'in' or 'mm'
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param ... parameters for the \code{\link[corrplot]{corrplot.mixed}} or
#' \code{\link[pheatmap]{pheatmap}} visualisation, for details see examples and
#'   help to corresponding functions
#'
#' @export
#'
#' @seealso \code{\link[pheatmap]{pheatmap}}, \code{\link[corrplot]{corrplot.mixed}}
plot_corr_matrix <- function(corr_matrix, flavor = 'corrplot', filename = NULL,
                             width = NA, height = NA, unit = c('cm','in','mm'),
                             plot_title = NULL, ...) {
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
#' @param peptide_annotation df with peptides and their corresponding proteins
#' @param protein_col the column name in \code{peptide_annotation} with protein names
#' @param peptide_col_name the column name in \code{peptide_annotation} with peptide names
#' @param flavor either corrplot from 'corrplot' package or heatmap, as in 'pheatmap'
#' @param filename path where the results are saved. If null the object is returned to the active window;
#' otherwise, the object is save into the file. Currently only pdf and png format is supported
#' @param width option  determining the output image width
#' @param height option  determining the output image width
#' @param unit units: 'cm', 'in' or 'mm'
#' @param plot_title The title of the plot
#' @param ... parameters for the corrplot visualisation
#'
#'
#' @export
#' @examples \donotrun{plot_corr_plot(q_norm_proteome, protein_name = 'Haao',
#'                peptide_annotation = peptide_annotation, prot.column = 'Gene',
#'                title = 'Haao protein peptides after quantile norm',
#'                number.cex=0.75, tl.cex = .75
#'                mar=c(0,0,1,0))}

#' @examples \donotrun{lower = "ellipse", upper = "number",
#'  tl.col = "black", diag = 'l', tl.pos = "lt", number.cex=0.75, tl.cex = .75}
#'
plot_protein_corrplot <- function(data_matrix,
           protein_name,
           peptide_annotation,
           protein_col = 'ProteinName',
           feature_id_col = 'peptide_group_label',
           flavor = 'corrplot',
           filename = NULL,
           width = NA, height = NA, unit = c('cm','in','mm'),
           plot_title = sprintf('Peptide correlation matrix of %s protein', protein_name), ...) {

    #extract peptides of the protein
    peptides = peptide_annotation %>%
      filter(UQ(sym(peptide_col_name)) %in% rownames(data_matrix)) %>%
      filter(UQ(sym(protein_col)) == protein_name) %>%
      pull(peptide_col_name) %>% as.character()
    data_matrix_sub = data_matrix[peptides,]
    corr_matrix = cor(t(data_matrix_sub), use = 'complete.obs')
    plot_corr_matrix(corr_matrix, plot_title = plot_title, flavor = flavor,
                     filename = filename, width = width, height = height, unit = unit, ...)
  }

#' Sample correlation matrix (heatmap)
#'
#' Plot correlation of selected samples
#'
#' @inheritParams proBatch
#' @param samples_to_plot string vector of samples in \code{data_matrix} to be used in the plot
#' @param filename path where the results are saved. If null the object is returned to the active window;
#' otherwise, the object is save into the file. Currently only pdf and png format is supported
#' @param width option  determining the output image width
#' @param height option  determining the output image width
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param ... parameters for the \code{\link[corrplot]{corrplot.mixed}} or
#' \code{\link[pheatmap]{pheatmap}} visualisation, for details see examples and
#'   help to corresponding functions
#'
#' @export
#'
#' @examples
#' ### Example 1: Plot heatmap of pre-specified samples
#' #'\dontrun{specified_samples = sample_annotation %>%
#' filter(RunID %in% paste('Run', 110:115, sep = '')) %>%
#' pull(FullRunName)
#' plot_samples_corr_heatmap(data_matrix, sample_to_plot = specified samples,
#'  flavor = 'pheatmap', cluster_rows = F, cluster_cols = F)
#'
#' }
#' ### Example 2: Plot corrplot of pre-specified samples
#' #' #'\dontrun{specified_samples = sample_annotation %>%
#' filter(RunID %in% paste('Run', 110:115, sep = '')) %>%
#' pull(FullRunName)
#' plot_samples_corr_heatmap(data_matrix, sample_to_plot = specified samples,
#'  flavor = 'corrplot', lower = "ellipse", upper = "number",
#'  tl.col = "black", diag = 'l', tl.pos = "lt", number.cex=0.75, tl.cex = .75)
#'
#' }
#' @seealso \code{\link[pheatmap]{pheatmap}}, \code{\link[corrplot]{corrplot.mixed}}
plot_sample_heatmap <- function(data_matrix, samples_to_plot = NULL,
                                      flavor = 'corrplot', filename = NULL,
                                      width = NA, height = NA, unit = c('cm','in','mm'),
                                      plot_title = sprintf('Correlation matrix of sample %s', samples_to_plot), ...){
  if(!is.null(samples_to_plot)){
    corr_matrix = cor(data_matrix[,samples_to_plot], use = 'complete.obs')
  } else {
    corr_matrix = cor(data_matrix, use = 'complete.obs')
  }
  plot_corr_matrix(corr_matrix, plot_title = plot_title, flavor = flavor,
                   filename = filename, width = width, height = height, unit = unit, ...)
}


#' Calculates correlation distribution for all pairs of the replicated samples
#'
#' @inheritParams proBatch
#' @param cor_proteome sample correlation matrix (square)
#' @param biospecimen_id_col column in `sample_annotation` that captures the
#' biological sample, that (possibly) was profiled several times as technical replicates.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column with code, such as the following:
#'  \code{sample_annotation %>% mutate(biospecimen_id = paste(tumorNormal, patientID))}
#' @param batch_col column in `sample_annotation` that should be used for
#'   batch comparison
#'
#' @return dataframe with the following columns, that are suggested to use for
#' plotting in \code{\link{plot_sample_corr_distribution}} as \code{plot_param}:
#' \enumerate{
#' \item \code{replicate}
#' \item \code{batch_the_same}
#' \item \code{batch_replicate}
#' \item \code{batches}
#' }
#' other columns are: \enumerate{
#' \item \code{sample_id_1} & \code{sample_id_2}, both generated from \code{sample_id_col} variable
#' \item \code{correlation} - correlation of two corresponding samples
#' \item \code{batch_1} & \code{batch_2} or analogous, created the same as \code{sample_id_1}
#' }
#'
#'
#' @export
#'
get_sample_corr_distrib <- function(cor_proteome, sample_annotation,
                                   sample_id_col = 'FullRunName',
                                   biospecimen_id_col = 'EarTag',
                                   batch_col = 'batch'){
  #since we only need unique pairs of samples, we create df of combinations to keep
  comb_to_keep = data.frame(t(combn(colnames(cor_proteome), 2)))
  names(comb_to_keep) = paste(sample_id_col, 1:2, sep = '_')

  spec_cols = c(biospecimen_id_col, batch_col)

  #transforming square matrix to long format
  corr_distribution = melt(cor_proteome,
                           varnames = paste(sample_id_col,1:2, sep = '_'),
                           value.name = 'correlation') %>%
    merge(comb_to_keep) %>%
    #merging with sample annotation, where we keep only sample_id_col, batch_col and biospecimen_col
    merge(sample_annotation %>% select(one_of(c(sample_id_col, spec_cols))),
          by.x = paste(sample_id_col,'1', sep = '_'),
          by.y = sample_id_col, all.x = T) %>%
    #to make it unambiguous, we rename columns, related to info of the left hand sample with "1" suffix, e.g. "Batch_1" meaning "batch of sample 1"
    data.table::setnames(old = spec_cols, new = paste(spec_cols, 1, sep = '')) %>%
    merge(sample_annotation %>% select(one_of(c(sample_id_col, spec_cols))),
          by.x = paste(sample_id_col,'2', sep = '_'),
          by.y = sample_id_col, all.x = T) %>%
    data.table::setnames(old = spec_cols, new = paste(spec_cols, 2, sep = '')) %>%
    #if biospecimen_1 and biospecimen_2 are the same, these samples are replicates
    mutate(replicate = (!!sym(paste(biospecimen_id_col,'1', sep = '')) ==
                          !!sym(paste(biospecimen_id_col,'2', sep = '')))) %>%
    #if batch_1 and batch_2 are the same, samples come from the same batch
    #some batches are more similar than the others, thus using 'batches' can reveal, that 'batch_1:batch_2' are more similar than 'batch_2:batch_3'
    mutate(batch_the_same = (!!sym(paste(batch_col,'1', sep = '')) ==
                               !!sym(paste(batch_col,'2', sep = ''))),
           batches = paste(!!sym(paste(batch_col,'1', sep = '')),
                           !!sym(paste(batch_col,'2', sep = '')), sep = ':')) %>%
    #to illustrate batch vs replicate interaction, this column captures all options
    mutate(batch_replicate = ifelse(replicate,
                                    ifelse(batch_the_same, 'same_batch\nsame_biospecimen', 'same_biospecimen\ndiff_batch'),
                                    ifelse(batch_the_same, 'same_batch\ndiff_biospecimen','diff_batch\ndiff_biospecimen')))
  return(corr_distribution)
}

#' Create violin plot of correlation distribution
#'
#' Useful to visualize within batch vs within replicate vs non-related sample correlation
#'
#' @param repeated_samples if `NULL`, only repeated sample correlation is plotted
#' @param biospecimen_id_col column in `sample_annotation` that captures the
#' biological sample, that (possibly) was profiled several times as technical replicates.
#'  Tip: if such ID is absent, but can be defined from several columns,
#'  create new \code{biospecimen_id} column with code, such as the following:
#'  \code{sample_annotation %>% mutate(biospecimen_id = paste(tumorNormal, patientID))}
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
#' @return \code{ggplot} type object with violin plot for each \code{plot_param}
#'
#' @export
#'
#' @examples
#' @seealso \code{\link{get_sample_corr_distrib}}, \code{\link[ggplot2]{ggplot}}
plot_sample_corr_distribution <- function(data_matrix, sample_annotation,
                                   repeated_samples = NULL,
                                   sample_id_col = 'FullRunName',
                                   batch_col = 'batch',
                                   biospecimen_id_col = 'EarTag',
                                   plot_title = 'Correlation distribution',
                                   plot_param = 'batch_replicate'){
  
  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(colnames(data_matrix))) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  
  corr_distribution <- function(data_matrix, repeated_samples, sample_annotation,
                                biospecimen_id_col, sample_id_col, batch_col) {
    if (!is.null(repeated_samples)){
      print('plotting correlation of repeated samples only')
      corr_matrix = cor(data_matrix[,repeated_samples], use = 'complete.obs')
    } else {
      corr_matrix = cor(data_matrix, use = 'complete.obs')
    }
    corr_distribution = get_sample_corr_distrib(cor_proteome = corr_matrix,
                                                sample_annotation = sample_annotation,
                                                sample_id_col = sample_id_col,
                                                biospecimen_id_col = biospecimen_id_col,
                                                batch_col = batch_col)
    return(corr_distribution)
  }
  if (!is.list(data_matrix)){
    corr_distribution = corr_distribution(data_matrix, repeated_samples, sample_annotation,
                                          biospecimen_id_col, sample_id_col, batch_col)
  } else {
    corr_distribution = lapply(1:length(data_matrix), function(i) {
      dm = data_matrix[[i]]
      corr_distribution = corr_distribution(dm, repeated_samples, sample_annotation,
                                            biospecimen_id_col, sample_id_col, batch_col)
      corr_distribution$Step = names(data_matrix)[i]
      return(corr_distribution)
    })
    corr_distribution = do.call(rbind, corr_distribution)
  }
  p <- ggplot(corr_distribution, aes_string(x = plot_param, y = 'correlation'))+
    geom_violin()+
    geom_boxplot(width = .1)+
    ggtitle(plot_title)+
    theme_classic()
  if (plot_param =='batches'){
    p = p + theme(axis.text.x = element_text(angle = 90))
  }
  if(is.list(data_matrix)){
    if(length(data_matrix) <= 4){
      p = p + theme_bw() + facet_grid(.~Step)
    }
    else {
      p = p + theme_bw()+facet_grid(Step ~ .)
    }
  }
  p = p + theme(plot.title = element_text(hjust = .5, face = 'bold'))
  return(p)
}




#' Transform square correlation matrix into long data frame of correlations
#'
#' @param peptide_cor peptide correlation matrix (square)
#' @param peptide_annotation df with petpides and their corresponding proteins
#' @param protein_col the column name in \code{peptide_annotation} with protein names
#' @param feature_id_col column in \code{peptide_annotation} that captures peptide names are found 
#'
#' @return
#' @export
#'
#' @examples
get_peptide_corr_df <- function(peptide_cor, peptide_annotation, protein_col = 'ProteinName',
                             feature_id_col = 'peptide_group_label'){
  #since we only need unique pairs of samples, we create df of combinations to keep
  comb_to_keep = data.frame(t(combn(colnames(peptide_cor), 2)))
  names(comb_to_keep) = paste(feature_id_col, 1:2, sep = '_')

  #transforming square matrix to long format
  corr_distribution = melt(peptide_cor,
                           varnames = paste(feature_id_col,1:2, sep = '_'),
                           value.name = 'correlation') %>%
    merge(comb_to_keep) %>%
    #merging with sample annotation, where we keep only sample_id_col, batch_col and biospecimen_col
    merge(peptide_annotation %>% select(one_of(c(feature_id_col, protein_col))),
          by.x = paste(feature_id_col,'1', sep = '_'),
          by.y = feature_id_col, all.x = T) %>%
    #to make it unambiguous, we rename columns, related to info of the left hand sample with "1" suffix, e.g. "Batch_1" meaning "batch of sample 1"
    data.table::setnames(old = protein_col, new = paste(protein_col, 1, sep = '')) %>%
    merge(peptide_annotation %>% select(one_of(c(feature_id_col, protein_col))),
          by.x = paste(feature_id_col,'2', sep = '_'),
          by.y = feature_id_col, all.x = T) %>%
    data.table::setnames(old = protein_col, new = paste(protein_col, 2, sep = '')) %>%
    #if protein_1 and protein_2 are the same
    mutate(same_protein = (!!sym(paste(protein_col,'1', sep = '')) ==
                          !!sym(paste(protein_col,'2', sep = ''))))
  return(corr_distribution)
}

#' Plot distribution of peptide correlations within one protein and between proteins
#'
#' @param data_matrix
#' @param peptide_annotation
#' @param protein_col
#' @param feature_id_col
#' @param plot_title
#' @param theme
#'
#' @return
#' @export
#'
#' @examples
plot_peptide_corr_distribution <- function(data_matrix, peptide_annotation,
                                        protein_col = 'ProteinName',
                                        feature_id_col = 'peptide_group_label',
                                        plot_title = 'Distribution of peptide correlation',
                                        theme = 'classic'){
  corr_distribution_prot <- function(data_matrix, peptide_annotation, protein_col, feature_id_col)
 {
    corr_matrix = cor(t(data_matrix), use = "pairwise.complete.obs")

    corr_distribution = get_peptide_corr_df(peptide_cor = corr_matrix,
                                            peptide_annotation = peptide_annotation,
                                            protein_col = protein_col,
                                            feature_id_col = feature_id_col)
    return(corr_distribution)
  }
  if (!is.list(data_matrix)){
    corr_distribution = corr_distribution_prot(data_matrix, peptide_annotation,
                                               protein_col, feature_id_col)
  } else {
    corr_distribution = lapply(1:length(data_matrix), function(i) {
      dm = data_matrix[[i]]
      corr_distribution = corr_distribution_prot(dm, peptide_annotation,
                                                 protein_col, feature_id_col)
      corr_distribution$Step = names(data_matrix)[i]
      return(corr_distribution)
    })
    corr_distribution = do.call(rbind, corr_distribution)%>%
      mutate(Step = factor(Step, levels = names(data_matrix)))
  }
  p <- ggplot(corr_distribution, aes_string(x = 'same_protein', y = 'correlation'))+
    geom_violin()+
    geom_boxplot(width = .1)+
    ggtitle(plot_title)+
    theme_classic()

  if(is.list(data_matrix)){
    if(length(data_matrix) <= 4){
      p = p + theme_bw() + facet_grid(.~Step)
    }
    else {
      p = p + theme_bw()+facet_grid(Step ~ .)
    }
  }
  p = p + theme(plot.title = element_text(hjust = .5, face = 'bold'))
  return(p)
}

#' Plot distribution of median correlations of peptides within same protein
#'
#' @param data_matrix
#' @param peptide_annotation
#' @param protein_col
#' @param feature_id_col
#' @param plot_title
#' @param theme
#'
#' @return
#' @export
#'
#' @examples
#'
#'
plot_within_prot_corr_distribution <- function(data_matrix_list, peptide_annotation,
                                               protein_col = 'ProteinName',
                                               feature_id_col = 'peptide_group_label',
                                               plot_title = 'Distribution of peptide correlation',
                                               theme = 'classic'){
  summarize_prot_corr_distribution <- function(data_matrix, peptide_annotation,
                                               protein_col, feature_id_col){
    corr_matrix = cor(t(data_matrix), use = "pairwise.complete.obs")

    corr_distribution = get_peptide_corr_df(peptide_cor = corr_matrix,
                                            peptide_annotation = peptide_annotation,
                                            protein_col = protein_col,
                                            feature_id_col = feature_id_col)
    corr_distribution = corr_distribution %>%
      filter(same_protein == T) %>%
      group_by(!!sym(paste(feature_id_col,'1', sep = '_'))) %>%
      summarize(median_correlation = median(correlation, na.rm = T))
    return(corr_distribution)
  }
  corr_distribution_list = lapply(1:length(data_matrix_list), function(i) {
    dm = data_matrix_list[[i]]
    corr_median_distribution = summarize_prot_corr_distribution(dm, peptide_annotation,
                                                         protein_col, feature_id_col)
    corr_median_distribution$Step = names(data_matrix_list)[i]
    return(corr_median_distribution)
  })
  corr_median_distribution_final = do.call(rbind, corr_distribution_list)
  corr_median_distribution_final = corr_median_distribution_final %>%
    mutate(Step = factor(Step, levels = names(data_matrix_list)))

  gg <- ggplot(corr_median_distribution_final, aes(y = median_correlation,
                                                   x = Step))+
    geom_violin()+
    geom_boxplot(width = .1)

  if(!is.null(plot_title))
    gg = gg + ggtitle(plot_title)
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }

  return(gg)
}

#' Plot distribution curves for each step to see if they shift by data normalization
#'
#' @param data_matrix_list
#' @param protein_name name of the protein, as specified in \code{protein_col}
#' of \code{peptide_annotation}
#' @param peptide_annotation
#' @param protein_col
#' @param feature_id_col
#' @param plot_title
#' @param theme
#'
#' @return
#' @export
#'
#' @examples
plot_peptide_correlation_distr_one_protein <- function(data_matrix_list,
                                                       protein_name,
                                                       peptide_annotation,
                                                       protein_col = 'ProteinName',
                                                       feature_id_col = 'peptide_group_label',
                                                       plot_title = sprintf('Distribution of peptide correlation at different correction steps,
                                                                            %s protein',protein_name),
                                                       theme = 'classic'){
  for(i in 2:length(data_matrix_list)){
    if (!all(rownames(data_matrix_list)[[1]] == rownames(data_matrix_list[[i]])))
        stop("rownames of different data matrices don't match, please check!")
  }
  feature_ids = rownames(data_matrix_list[[1]])
  peptides = peptide_annotation %>%
    filter(UQ(sym(feature_id_col)) %in% feature_ids) %>%
    filter(UQ(sym(protein_col)) == protein_name) %>%
    pull(feature_id_col) %>% as.character()

  data_matrix_list = lapply(data_matrix_list, function(dm) {
    dm_sub = dm[peptides,]
    return(dm_sub)
  })

  corr_distribution_prot <- function(data_matrix, peptide_annotation, protein_col, feature_id_col)
  {
    corr_matrix = cor(t(data_matrix), use = 'complete.obs')

    corr_distribution = get_peptide_corr_df(peptide_cor = corr_matrix,
                                            peptide_annotation = peptide_annotation,
                                            protein_col = protein_col,
                                            feature_id_col = feature_id_col)
    return(corr_distribution)
  }
  corr_distribution = lapply(1:length(data_matrix_list), function(i) {
    dm = data_matrix_list[[i]]
    corr_distribution = corr_distribution_prot(dm, peptide_annotation, protein_col,
                                          feature_id_col)
    corr_distribution$Step = names(data_matrix_list)[i]
    return(corr_distribution)
  })

  corr_distribution = do.call(rbind, corr_distribution)
  gg <- ggplot(corr_distribution, aes(linetype = Step, x = correlation))+
    geom_line(stat="density")
  if(!is.null(plot_title))
    gg = gg + ggtitle(plot_title)
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  return(gg)
}
