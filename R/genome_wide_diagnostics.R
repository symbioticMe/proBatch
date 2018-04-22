#' Plot the sample average
#'
#' @param data_matrix
#' @param sample_annotation
#' @param sample_id_col
#' @param batch_column
#' @param order_column
#' @param color_by_batch
#' @param theme
#' @param title
#' @param color_scheme
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
plot_sample_mean <- function(data_matrix, sample_annotation,
                             sample_id_col = 'FullRunName',
                             order_column = 'order',
                             batch_column = NULL,
                             color_by_batch = F, theme = 'classic',
                             title = NULL, color_scheme = 'brewer'){
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
#'
#' @param list_of_dfs
#' @param sample_annotation
#' @param batch_column
#' @param steps
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
boxplot_all_steps <- function(list_of_dfs, sample_annotation, batch_column,
                              steps = NULL){
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
