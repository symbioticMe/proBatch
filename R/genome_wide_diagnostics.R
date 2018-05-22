#' Plot sample average or distribution (boxplot) vs order (if the real running order available)
#' for quick assessment of trends associated, overall or specific covariate-associated (see `batch_column` and `facet_column`)
#' @param data_matrix features (in rows) vs samples (in columns) matrix,
#' with feature IDs in rownames and file/sample names as colnames.
#' in most function, it is assumed that this is the log transformed version of the original data
#' @param df_long data frame where each row is a single feature in a single sample,
#' thus it has minimally, `sample_id_col`, `feature_id_column` and `measure_column`,
#' but usually also `m_score` (in OpenSWATH output result file)
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names) 2) biological and 3) technical covariates (batches etc)
#' @param sample_id_col name of the column in sample_annotation file,
#' where the filenames (colnames of the data matrix are found)
#' @param measure_column if `df_long` is among the parameters, it is the column with expression/abundance/intensity,
#' otherwise, it is used internally for consistency
#' @param batch_column column in `sample_annotation` that should be used for batch comparison
#' @param order_column column where running order is specified.
#' @param color_by_batch should the each batch be represented with its own color?
#' @param color_scheme named vector, names corresponding to unique batch values as specified in `sample_annotation`
#' @param facet_column recommended if more than one batch covariate is present. Faceting is most suited to examine instruments separately
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see examples)
#' @param title Title of the plot (usually, processing step + representation level (fragments, transitions, proteins))
#' @return ggplot2 class object. Thus, all aesthetics can be overriden
#' @name plot_sample_means_or_boxplots


#' @name plot_sample_means_or_boxplots
#'
#' @export
#' @import ggplot2
#' @import dplyr
#' @import rlang
#'
#' @examples
plot_sample_mean <- function(data_matrix, sample_annotation = NULL,
                             sample_id_col = 'FullRunName',
                             order_column = 'order',
                             batch_column = NULL,
                             facet_column = 'instrument',
                             color_by_batch = F, color_scheme = 'brewer',
                             theme = 'classic',
                             title = NULL){
  sample_average = colMeans(data_matrix)
  names(sample_average) = colnames(data_matrix)

  df_ave = data.frame(average = sample_average,
                      order_temp_col = 1:length(sample_average),
                      sample_id_col = colnames(data_matrix))
  names(df_ave)[names(df_ave) == "sample_id_col"] <- sample_id_col
  df_ave = df_ave %>%
    merge(sample_annotation, by = sample_id_col)
  if (!(order_column %in% names(sample_annotation))){
    warning('order column not found in sample annotation, taking order of files in the data matrix instead')
    order_column = 'order_temp_col'
  }
  gg = ggplot(df_ave, aes_string(x = order_column, y = 'average'))+
                geom_point()
  if(color_by_batch & !is.null(batch_column)){
    gg = gg + aes_string(color = batch_column)
  }
  if(!is.null(batch_column)){
    #TODO: modify this for multi-instrument case
    batch.tipping.points = cumsum(table(sample_annotation[[batch_column]]))+.5
    gg = gg + geom_vline(xintercept = batch.tipping.points,
                         color = 'grey', linetype = 'dashed')
  }
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if(!is.null(title)) gg = gg + ggtitle(title)+
    theme(plot.title = element_text(face = 'bold',hjust = .5))
  if(length(color_scheme) == 1 & color_scheme == 'brewer'){
    gg = gg + scale_color_brewer(palette = 'Set1')
  } else{
    gg = gg + scale_color_manual(values = color_scheme)
  }

  return(gg)
}

#' plot boxplot of data, optionally colored by batch
#'
#' @param data_df_long
#' @param sample_annotation
#' @param batch_column
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
gg_boxplot <- function(data_df_long, sample_annotation, batch_column,
                      order_column = 'order', measure_col = 'Intensity',
                      fill_batch = T, theme = 'classic', title = NULL){
  if (!all(names(sample_annotation) %in% names(data_df_long))){
    data_df_long = data_df_long %>% merge(sample_annotation)
  }
  gg = ggplot(data_df_long, aes_string(x = order_column, y = measure_col,
                                       group = order_column))+
    geom_boxplot(outlier.size = .15)+
    theme_bw()+
    theme(plot.title = element_text(face = 'bold', hjust = .5),
          legend.position="top", axis.text.x = element_text(angle = 90))
  if(fill_batch){
    gg = gg + aes_string(fill = batch_column)
  }
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if (!is.null(title)){
    gg = gg + ggtitle(title)+theme(plot.title = element_text(hjust = .5, face = 'bold', size = 16))
  }
  return(gg)
}

#' Plot boxplots to compare various data normalization steps/approaches
#' WARNING: extremely slow for big dataframes
#'
#' @param list_of_dfs list of data frames of format, specified in `plot_boxplot`
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names) 2) biological and 3) technical covariates (batches etc)
#' @param batch_column column in `sample_annotation` that should be used for batch comparison
#' @param step normalization step (e.g. `Raw` or `Quantile_normalized` or `qNorm_ComBat`).
#' Useful if consecutive steps are compared in plots.
#' Note that in plots these are usually ordered alphabetically, so it's worth naming with numbers, e.g. `1_raw`, `2_quantile`

#'
#' @return ggplot object
#' @export
#' @import tidyverse
#'
#' @examples
#' @seealso \code{\link{plot_boxplot}}
boxplot_all_steps <- function(list_of_dfs, sample_annotation, batch_column,
                              step = NULL){
  if(`*`(dim(list_of_dfs[[1]])[1], dim(list_of_dfs[[1]])[2]) * length(list_of_dfs) > 5*10^5){
    warning('Data matrices are huge, be patient, this might take a while (or crash)')
  }
  add_processing_step <- function(i, list_of_dfs, steps) {
    df = list_of_dfs[[i]]
    df$step = steps[i]
    list_of_dfs[[i]] = df
  }
  if (!is.null(steps) |
      (any(sapply(list_of_dfs, function(df) (!'step' %in% names(df)))))){
    list_of_dfs = lapply(1:length(list_of_dfs), add_processing_step, list_of_dfs, steps)
  }

  joined_proteome = do.call(rbind, list_of_dfs)
  gg = gg_boxplot(joined_proteome, sample_annotation, batch_column) +
    facet_grid(.~step)
  return(gg)
}

#' cluster the data matrix to visually inspect which confounder dominates
#'
#' @param proteome
#' @param color_df
#' @param title
#'
#' @return
#' @export
#' @importFrom WGCNA plotDendroAndColors
#'
#' @examples
cluster_samples <- function(data_matrix, color_df, plot_title, ...){
  dist_matrix = dist(t(as.matrix(data_matrix)))
  hierarchical_clust = hclust(dist_matrix)
  plotDendroAndColors(hierarchical_clust, color_df, rowTextAlignment = 'left',
                      main = plot_title,
                      hang = -0.1, addGuide = T, ...)
}

#' Quantify variance distribution by variable
#'
#' @param data_matrix
#' @param sample_annotation
#' @param technical_covariates
#' @param biological_covariates
#' @param plot_title
#' @param colors_for_bars
#' @param threshold_pca
#' @param threshold_var
#'
#' @return
#' @export
#' @import ggplot2
#' @importFrom pvca pvcaBatchAssess
#'
#' @examples
plot_pvca <- function(data_matrix, sample_annotation, sample_id_column = 'FullRunName',
                 technical_covariates = NULL, biological_covariates = NULL,
                 plot_title = NULL, colors_for_bars = NULL, threshold_pca = .6,
                 threshold_var = .01, theme = 'classic'){
  factors_for_PVCA = c(technical_covariates, biological_covariates)
  if (!is.null(sample_id_column)){
    if(sample_id_column %in% names(sample_annotation)){
      rownames(sample_annotation) = sample_annotation[[sample_id_column]]
    }
  }

  tech_interactions = expand.grid(technical_covariates, technical_covariates) %>%
    mutate(tech_interaction = paste(Var1, Var2, sep = ':')) %>%
    pull(tech_interaction)
  biol_interactions = expand.grid(biological_covariates, biological_covariates) %>%
    mutate(biol_interactions = paste(Var1, Var2, sep = ':')) %>%
    pull(biol_interactions)

  if(!(all(rownames(sample_annotation) %in% colnames(data_matrix)))){
    stop('check data matrix column names or these in sample annotation')
  }
  covrts.annodf = Biobase::AnnotatedDataFrame(data=sample_annotation)
  expr_set = Biobase::ExpressionSet(data_matrix[,rownames(sample_annotation)], covrts.annodf)
  pvcaAssess <- pvcaBatchAssess (expr_set, factors_for_PVCA, threshold = threshold_pca)
  pvcaAssess_df = data.frame(weights = as.vector(pvcaAssess$dat),
                             label = pvcaAssess$label,
                             stringsAsFactors = F)

  label_of_small = sprintf('Below %1.0f%%', 100*threshold_var)
  if (sum(pvcaAssess_df$weights < threshold_var) > 1){
    pvca_res_small = sum(pvcaAssess_df$weights[pvcaAssess_df$weights < threshold_var])
    pvca_res = pvcaAssess_df[pvcaAssess_df$weights >= threshold_var, ]
    pvca_res_add = data.frame(weights = pvca_res_small, label = label_of_small)
    pvca_res = rbind(pvca_res, pvca_res_add)
  }
  else {pvca_res = pvcaAssess_df}

  technical_covariates = c(technical_covariates, tech_interactions)
  biological_covariates = c(biological_covariates, biol_interactions)
  pvca_res = pvca_res %>% mutate(category = ifelse(label %in% technical_covariates, 'technical',
                                                   ifelse(label %in% biological_covariates, 'biological',
                                                          ifelse(label %in% c(label_of_small, 'resid'), 'residual', 'biol:techn'))))

  pvca_res = pvca_res %>%
    arrange(desc(weights)) %>%
    arrange(label == label_of_small) %>%
    arrange(label == 'resid')

  pvca_res$label = factor(pvca_res$label, levels = pvca_res$label)

  y_title = 'Weighted average proportion variance'
  gg  = ggplot(pvca_res, aes(x = label, y = weights, fill = category))+
    geom_bar(stat = 'identity', color = 'black', size = 1.5)+
    ylab(y_title)
  if(theme == 'classic'){
    gg = gg +theme_classic()
  }
  if(is.null(colors_for_bars)){
    colors_for_bars = c('grey', wesanderson::wes_palettes$Rushmore[3:5])
    names(colors_for_bars) = c('residual', 'biological', 'biol:techn', 'technical')
    gg = gg + scale_fill_manual(values = colors_for_bars)
  }
  gg = gg +
    theme(axis.title.x = NULL, axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
    xlab(NULL)+
    theme(text = element_text(size=15))+
    guides(fill=guide_legend(override.aes=list(color=NA), title=NULL))
  return(list(plot =gg, df = pvca_res))
}

#' plot PCA plot
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix,
#' with feature IDs in rownames and file/sample names as colnames.
#' in most function, it is assumed that this is the log transformed version of the original data
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names) 2) biological and 3) technical covariates (batches etc)
#' @param color_by column name (as in `sample_annotation`) to color by
#' @param PC_to_plot principal component numbers for x and y axis
#' @param colors_for_factor named vector of colors for the `color_by` variable
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see examples)
#'
#' @return ggplot scatterplot colored by factor levels of column specified in `factor_to_color`
#' @export
#' @import ggplot2
#' @import ggfortify
#'
#' @examples
plot_pca <- function(data_matrix, sample_annotation, factor_to_color,
                     PC_to_plot = c(1,2), colors_for_factor = NULL, theme = 'classic'){
  gg = autoplot(prcomp(t(data_matrix)), data = sample_annotation, colour = factor_to_color,
           x = PC_to_plot[1], y = PC_to_plot[2])
  if (theme == 'classic'){
    gg = gg + theme_classic()
  }
  if(!is.null(colors_for_factor)){
    gg = gg + scale_color_manual(values = colors_for_factor)
  }
  return(gg)
}

#' Plot the heatmap of samples
#'
#' @param data_matrix
#' @param sample_annotation
#' @param fill_the_missing
#' @param cluser_rows
#' @param cluster_cols
#' @param annotation_color_list
#' @param filename
#' @param ...
#'
#' @return
#' @export
#' @import pheatmap
#'
#' @examples
plot_heatmap <- function(data_matrix, sample_annotation, fill_the_missing = T,
                         cluser_rows = F, cluster_cols = F,
                         annotation_color_list = NA,
                         filename = NA,
                         ...){
  if(fill_the_missing) {
    data_matrix[is.na(data_matrix)] = 0
    heatmap_color = c('black', colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100))
  }
  else {
    heatmap_color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  }
  p <- pheatmap(data_matrix, cluster_rows = cluser_rows, cluster_cols = cluster_cols,
           color = heatmap_color,
           annotation_col = sample_annotation, annotation_colors = color_list,
           filename = filename, ...)
  return(p)
}
