#' Plot per-sample average or boxplot (distribution) vs order (if the real
#' running order available)
#' @details functions for quick visual assessment of trends associated, overall
#'   or specific covariate-associated (see `batch_column` and `facet_column`)
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param df_long data frame where each row is a single feature in a single
#'   sample, thus it has minimally, `sample_id_col`, `feature_id_column` and
#'   `measure_column`, but usually also `m_score` (in OpenSWATH output result
#'   file)
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_column if `df_long` is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @param batch_column column in `sample_annotation` that should be used for
#'   batch comparison
#' @param order_column column where running order is specified.
#' @param color_by_batch should the each batch be represented with its own
#'   color?
#' @param color_scheme named vector, names corresponding to unique batch values
#'   as specified in `sample_annotation`
#' @param facet_column recommended if more than one batch covariate is present.
#'   Faceting is most suited to examine instruments separately
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see
#'   examples)
#' @param title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param order_per_facet if order is defined ignoring facets (usually
#'   instrument), re-define order per-batch
#' @return ggplot2 class object. Thus, all aesthetics can be overriden
#' @name plot_sample_means_or_boxplots

#' @name plot_sample_means_or_boxplots
#'
#' @export
#'
#' @examples
plot_sample_mean <- function(data_matrix, sample_annotation = NULL,
                             sample_id_col = 'FullRunName',
                             order_column = 'order',
                             batch_column = NULL,
                             facet_column = 'instrument',
                             color_by_batch = F, color_scheme = 'brewer',
                             theme = 'classic',
                             title = NULL, order_per_facet = F){
  sample_average = colMeans(data_matrix, na.rm = T)
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
  if(!is.null(facet_column)){
    if(!(facet_column %in% names(df_ave))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present in the data,
                   check sample annotation data frame', facet_column))
    }
    if (order_per_facet){
      df_ave = df_ave %>%
        group_by_at(vars(one_of(facet_column))) %>%
        mutate(order = rank(UQ(sym(order_column))))
    }
  }
  gg = ggplot(df_ave, aes_string(x = order_column, y = 'average'))+
                geom_point()
  if(color_by_batch & !is.null(batch_column)){
    gg = gg + aes_string(color = batch_column)
    if(length(color_scheme) == 1 & color_scheme == 'brewer'){
      n_batches <- length(unique(sample_annotation[[batch_column]]))
      if (n_batches <= 9){
        gg = gg + scale_color_brewer(palette = 'Set1')
      } else {
        if (n_batches <= 12){
          gg = gg + scale_color_brewer(palette = 'Set3')
        } else {
          warning(sprintf('brewer palettes have maximally 12 colors, you specified %s batches,
                consider defining color scheme with sample_annotation_to_colors function', n_batches))
        }
      }

    } else{
      gg = gg + scale_color_manual(values = color_scheme)
    }
  }
  if(!is.null(batch_column)){
    if (!is.null(facet_column)){
      order_vars <- c(facet_column, order_column)
      batch_vars = c(facet_column, batch_column)
      tipping.points = df_ave %>%
        arrange(!!!syms(order_vars))%>%
        group_by(!!!syms(batch_vars)) %>%
        summarise(batch_size = n()) %>%
        group_by(!!sym(facet_column)) %>%
        mutate(tipping.points = cumsum(batch_size))%>%
        mutate(tipping.poings = tipping.points+.5)
      gg = gg + geom_vline(data = tipping.points, aes(xintercept = tipping.poings),
                           color = 'grey', linetype = 'dashed')
    } else {
      batch.tipping.points = cumsum(table(sample_annotation[[batch_column]]))+.5
      gg = gg + geom_vline(xintercept = batch.tipping.points,
                           color = 'grey', linetype = 'dashed')
    }
  }
  if(!is.null(facet_column)){
    gg = gg + facet_wrap(as.formula(paste("~", facet_column)),
                         dir = 'v', scales = "free_x")
  }

  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if(!is.null(title)) gg = gg + ggtitle(title)+
    theme(plot.title = element_text(face = 'bold',hjust = .5))


  return(gg)
}

#' @name plot_sample_means_or_boxplots
#'
#' @export
#'
#' @examples
plot_boxplot <- function(df_long, sample_annotation = NULL,
                       sample_id_column = 'FullRunName',
                       measure_col = 'Intensity',
                       order_column = 'order',
                       batch_column = 'MS_batch.final',
                       facet_column = 'instrument',
                       color_by_batch = T, color_scheme = 'brewer',
                       theme = 'classic',
                       title = NULL, order_per_facet = F){
  if (!all(c(batch_column, sample_id_column) %in% names(df_long))){
    if (!is.null(sample_annotation)){
      df_long = df_long %>% merge(sample_annotation,
                                            by = sample_id_column)
    } else {
      if (color_by_batch){
        stop('batches cannot be colored if the batch column cannot be defined,
             check sample_annotation and data matrix')
      }
    }
  }
  if (is.null(order_column)){
    warning('order column not defined, taking order of files in the data matrix instead')
    order_column = 'order_temp_col'
    df_long[[order_column]] = match(df_long[[sample_id_column]],
                                         unique(df_long[[sample_id_column]]))
  } else if (!(order_column %in% names(sample_annotation)) &
             !(order_column %in% names(df_long))){
    warning('order column not found in sample annotation, taking order of files in the data matrix instead')
    order_column = 'order_temp_col'
    df_long[[order_column]] = match(df_long[[sample_id_column]],
                                         unique(df_long[[sample_id_column]]))
    order_per_facet = T
  }

  if (order_per_facet){
    if (!is.null(facet_column)){
      warning('defining order within each facet')
      df_long = df_long %>%
        group_by_at(vars(one_of(facet_column))) %>%
        mutate(order = rank(UQ(sym(order_column))))
    }
  }

  gg = ggplot(df_long, aes_string(x = order_column, y = measure_col,
                                       group = order_column))+
    geom_boxplot(outlier.size = .15)+
    theme_bw()+
    theme(plot.title = element_text(face = 'bold', hjust = .5))


  if (!is.numeric(df_long[[order_column]])){
    warning(sprintf('order column is not numeric, assuming the order is irrelevant
                    of %s order follows the run order.', order_column))
    df_long[[order_column]] = factor(df_long[[order_column]],
                                              levels = df_long[[order_column]])
    gg = gg +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  }
  if(color_by_batch){
    gg = gg + aes_string(fill = batch_column)
    if(length(color_scheme) == 1 & color_scheme == 'brewer'){
      n_batches <- length(unique(df_long[[batch_column]]))
      if(n_batches < 9){
        gg = gg + scale_fill_brewer(palette = 'Set1')

      } else {
          if (n_batches <= 12){
            gg = gg + scale_fill_brewer(palette = 'Set3')
          } else {
            warning(sprintf('brewer palettes have maximally 12 colors, you specified %s batches,
                consider defining color scheme with sample_annotation_to_colors function', n_batches))
      }
      }
    } else{
      gg = gg + scale_fill_manual(values = color_scheme)
    }
  }
  if(!is.null(facet_column)){
    if(!(facet_column %in% names(df_long))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present in the data,
                   check sample annotation data frame', facet_column))
    }
    gg = gg + facet_wrap(as.formula(paste("~", facet_column)),
                         dir = 'v', scales = "free_x")
  }
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if (!is.null(title)){
    gg = gg + ggtitle(title)+
      theme(plot.title = element_text(hjust = .5, face = 'bold', size = 16))
  }
  if (max(df_long[[order_column]]) > 100){
    gg = gg + theme(legend.position="top")
  }
  return(gg)
}

#' Plot boxplots to compare various data normalization steps/approaches WARNING:
#' extremely slow for big dataframes
#'
#' @param list_of_dfs list of data frames of format, specified in `plot_boxplot`
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param batch_column column in `sample_annotation` that should be used for
#'   batch comparison
#' @param step normalization step (e.g. `Raw` or `Quantile_normalized` or
#'   `qNorm_ComBat`). Useful if consecutive steps are compared in plots. Note
#'   that in plots these are usually ordered alphabetically, so it's worth
#'   naming with numbers, e.g. `1_raw`, `2_quantile`

#'
#' @return ggplot object
#' @export
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
  if (!is.null(step) |
      (any(sapply(list_of_dfs, function(df) (!'step' %in% names(df)))))){
    list_of_dfs = lapply(1:length(list_of_dfs), add_processing_step, list_of_dfs, step)
  }

  joined_proteome = do.call(rbind, list_of_dfs)
  gg = gg_boxplot(joined_proteome, sample_annotation, batch_column) +
    facet_grid(.~step)
  return(gg)
}

#' cluster the data matrix to visually inspect which confounder dominates
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param distance distance metric used for clustering
#' @param agglomeration agglomeration methods as used by `hclust`
#' @param color_df data frame of colors, as created by
#'   `sample_annotation_to_colors`
#' @param ... other parameters of `plotDendroAndColors` from `WGCNA` package
#'
#' @export
#'
#' @examples
#' @seealso \code{\link{hclust}}, \code{\link{sample_annotation_to_colors}},
#'   \code{\link{plotDendroAndColors}}
plot_clustering <- function(data_matrix, color_df, title = 'Clustering of raw samples',
                            distance = "euclidean",
                            agglomeration = 'complete',
                            label_samples = T, label_font = .2,
                            plot_title = NULL,
                              ...){
  dist_matrix = dist(t(as.matrix(data_matrix)), method = distance)
  hierarchical_clust = hclust(dist_matrix, method = agglomeration)
  if (label_samples){
    if (ncol(data_matrix) > 80){
      warning('Too many samples, adjust the font with `label_font` argument or
              remove labels by setting `label_samples = F` in function call')
    }
    WGCNA::plotDendroAndColors(hierarchical_clust, color_df, rowTextAlignment = 'left',
                        main = plot_title,
                        hang = -0.1, addGuide = T, dendroLabels = NULL,
                        cex.dendroLabels = label_font, ...)

  } else{
    WGCNA::plotDendroAndColors(hierarchical_clust, color_df, rowTextAlignment = 'left',
                        main = plot_title,
                        hang = -0.1, addGuide = T, dendroLabels = F, ...)
  }

}

PVCA <- function(data_matrix, sample_annotation, factors_for_PVCA, threshold_pca, threshold_var = Inf) {
  covrts.annodf = Biobase::AnnotatedDataFrame(data=sample_annotation)
  expr_set = Biobase::ExpressionSet(data_matrix[,rownames(sample_annotation)], covrts.annodf)
  pvcaAssess <- pvca::pvcaBatchAssess (expr_set, factors_for_PVCA, threshold = threshold_pca)
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
  return(pvca_res)
}

#' Plot variance distribution by variable
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param technical_covariates vector `sample_annotation` column names that are
#'   technical covariates
#' @param biological_covariates vector `sample_annotation` column names, that
#'   are biologically meaningful covariates
#' @param title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#' @param threshold_pca the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param threshold_var the percentile value of weight each of the covariates
#'   needs to explain (the rest will be lumped together)
#'
#' @return list of two items: plot =gg, df = pvca_res
#' @export
#'
#' @examples
#' @seealso \code{\link{sample_annotation_to_colors}}
plot_pvca <- function(data_matrix, sample_annotation, sample_id_column = 'FullRunName',
                 technical_covariates = c('MS_batch', 'instrument'),
                 biological_covariates = c('cell_line','drug_dose'),
                 title = NULL, colors_for_bars = NULL, threshold_pca = .6,
                 threshold_var = .01, theme = 'classic'){
  factors_for_PVCA = c(technical_covariates, biological_covariates)
  if (!is.null(sample_id_column)){
    if(sample_id_column %in% names(sample_annotation)){
      rownames(sample_annotation) = sample_annotation[[sample_id_column]]
    }
  }

  pvca_res = PVCA(data_matrix, sample_annotation, factors_for_PVCA,
                   threshold_pca, threshold_var = threshold_var)

  tech_interactions = expand.grid(technical_covariates, technical_covariates) %>%
    mutate(tech_interactions = paste(Var1, Var2, sep = ':')) %>%
    pull(tech_interactions)
  biol_interactions = expand.grid(biological_covariates, biological_covariates) %>%
    mutate(biol_interactions = paste(Var1, Var2, sep = ':')) %>%
    pull(biol_interactions)

  if(!(all(rownames(sample_annotation) %in% colnames(data_matrix)))){
    stop('check data matrix column names or these in sample annotation')
  }

  technical_covariates = c(technical_covariates, tech_interactions)
  biological_covariates = c(biological_covariates, biol_interactions)
  pvca_res = pvca_res %>% mutate(category = ifelse(label %in% technical_covariates, 'technical',
                                                   ifelse(label %in% biological_covariates, 'biological',
                                                          ifelse(label %in% c(label_of_small, 'resid'), 'residual', 'biol:techn'))))

  #ToDo: create a ranking function for PVCA items, so that these are plotted by the following logic and not by alphabetical order
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

  }
  gg = gg + scale_fill_manual(values = colors_for_bars)

  if (!is.null(title)){
    gg = gg + ggtitle(title)
  }

  #cosmetic updates to the plot
  gg = gg +
    theme(axis.title.x = NULL, axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
    xlab(NULL)+
    theme(text = element_text(size=15))+
    guides(fill=guide_legend(override.aes=list(color=NA), title=NULL))

  return(list(plot =gg, df = pvca_res))
}

#' plot PCA plot
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param color_by column name (as in `sample_annotation`) to color by
#' @param PC_to_plot principal component numbers for x and y axis
#' @param colors_for_factor named vector of colors for the `color_by` variable
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see
#'   examples)
#'
#' @return ggplot scatterplot colored by factor levels of column specified in
#'   `factor_to_color`
#' @export
#'
#' @examples
plot_pca <- function(data_matrix, sample_annotation, color_by = 'MS_batch',
                     PC_to_plot = c(1,2),
                     colors_for_factor = NULL, theme = 'classic'){
  if(length(color_by) > 1){
    warning('Coloring by the first column specified')
    color_by = color_by[1]
  }
  pr_comp_res <- prcomp(t(data_matrix))
  gg = autoplot(pr_comp_res, data = sample_annotation,
                colour = color_by,
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
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param fill_the_missing boolean value determining if missing values should be
#'   substituted with -1 (and colored with black)
#' @param cluster_rows boolean value determining if rows should be clustered
#' @param cluster_cols boolean value determining if columns should be clustered
#' @param annotation_color_list list specifying colors for columns (samples).
#'   Best created by `sample_annotation_to_colors`
#' @param filename filepath where to save the image
#' @param ... other parameters of pheatmap
#'
#' @return object returned by `pheatmap`
#' @export
#'
#' @examples
#' @seealso \code{\link{sample_annotation_to_colors}}, \code{\link{pheatmap}}
plot_heatmap <- function(data_matrix, sample_annotation = NULL, fill_the_missing = T,
                         cluster_rows = F, cluster_cols = F,
                         annotation_color_list = NA,
                         filename = NA,
                         ...){
  if(fill_the_missing) {
    data_matrix[is.na(data_matrix)] = -1
    heatmap_color = c('black', colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100))
  }
  else {
    heatmap_color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100)
  }
  if (is.null(sample_annotation)){
    sample_annotation = NA
  }
  p <- pheatmap(data_matrix, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
           color = heatmap_color,
           annotation_col = sample_annotation, annotation_colors = annotation_color_list,
           filename = filename, ...)
  return(p)
}
