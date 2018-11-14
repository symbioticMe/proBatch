#' Plot per-sample mean or boxplot (showing median and quantiles) vs order (if the real
#' running order available)
#' @details functions for quick visual assessment of trends associated, overall
#'   or specific covariate-associated (see `batch_col` and `facet_col`)
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param df_long data frame where each row is a single feature in a single
#'   sample, thus it has minimally, `sample_id_col`, `feature_id_col` and
#'   `measure_col`, but usually also `m_score` (in OpenSWATH output result
#'   file)
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if `df_long` is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @param batch_col column in `sample_annotation` that should be used for
#'   batch comparison
#' @param order_col column where running order is specified.
#' @param color_by_batch should the each batch be represented with its own
#'   color?
#' @param color_scheme named vector, names corresponding to unique batch values
#'   as specified in `sample_annotation`
#' @param facet_col recommended if more than one batch covariate is present.
#'   Faceting is most suited to examine instruments separately
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see
#'   examples)
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param order_per_facet if order is defined ignoring facets (usually
#'   instrument), re-define order per-batch
#' @return ggplot2 class object. Thus, all aesthetics can be overriden
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @name plot_sample_means_or_boxplots

#' @name plot_sample_mean_or_boxplot
#'
#' @export
#'
#' @examples
plot_sample_mean <- function(data_matrix, sample_annotation = NULL,
                             sample_id_col = 'FullRunName',
                             order_col = 'order',
                             batch_col = NULL,
                             facet_col = 'instrument',
                             color_by_batch = F, color_scheme = 'brewer',
                             theme = 'classic',
                             plot_title = NULL, order_per_facet = F,
                             vline_color = 'grey',
                             ylimits = NULL){
  sample_average = colMeans(data_matrix, na.rm = T)
  names(sample_average) = colnames(data_matrix)
  
  df_ave = data.frame(Average_Intensity = sample_average,
                      order_temp_col = 1:length(sample_average),
                      sample_id_col = colnames(data_matrix))
  names(df_ave)[names(df_ave) == "sample_id_col"] <- sample_id_col

  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(df_ave[[sample_id_col]])) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  df_ave = df_ave %>%
    merge(sample_annotation, by = sample_id_col)
  
  sample_annotation = sample_annotation %>%
    subset(sample_annotation[[sample_id_col]] %in% df_ave[[sample_id_col]])
  
  if (!(order_col %in% names(sample_annotation))){
    warning('order column not found in sample annotation, taking order of files in the data matrix instead')
    order_col = 'order_temp_col'
  }
  if(!is.null(facet_col)){
    if(!(facet_col %in% names(df_ave))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present in the data,
                   check sample annotation data frame', facet_col))
    }
    if (order_per_facet){
      df_ave = df_ave %>%
        group_by_at(vars(one_of(facet_col))) %>%
        mutate(order = rank(UQ(sym(order_col))))
    }
  }
  gg = ggplot(df_ave, aes_string(x = order_col, y = 'Average_Intensity'))+
    geom_point()
  if(!is.null(ylimits)){
    gg = gg +
      ylim(ylimits)
  }
  
  if(color_by_batch & !is.null(batch_col)){
    gg = gg + aes_string(color = batch_col)
    if(length(color_scheme) == 1 & color_scheme == 'brewer'){
      n_batches <- length(unique(sample_annotation[[batch_col]]))
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
  if(!is.null(batch_col)){
    if (!is.null(facet_col)){
      order_vars <- c(facet_col, order_col)
      batch_vars = c(facet_col, batch_col)
      tipping.points = df_ave %>%
        arrange(!!!syms(order_vars))%>%
        group_by(!!!syms(batch_vars)) %>%
        summarise(batch_size = n()) %>%
        group_by(!!sym(facet_col)) %>%
        mutate(tipping.points = cumsum(batch_size))%>%
        mutate(tipping.poings = tipping.points+.5)
      gg = gg + geom_vline(data = tipping.points, aes(xintercept = tipping.poings),
                           color = vline_color, linetype = 'dashed')
    } else {
      batch.tipping.points = cumsum(table(sample_annotation[[batch_col]]))+.5
      gg = gg + geom_vline(xintercept = batch.tipping.points,
                           color = vline_color, linetype = 'dashed')
    }
  }
  if(!is.null(facet_col)){
    gg = gg + facet_wrap(as.formula(paste("~", facet_col)),
                         dir = 'v', scales = "free_x")
  }
  
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if(!is.null(plot_title)) gg = gg + ggtitle(plot_title)+
    theme(plot.title = element_text(face = 'bold',hjust = .5))
  
  
  return(gg)
}


#' @name plot_sample_mean_or_boxplot
#'
#' @export
#'
#' @examples
plot_boxplot <- function(df_long, sample_annotation = NULL,
                       sample_id_col = 'FullRunName',
                       measure_col = 'Intensity',
                       order_col = 'order',
                       batch_col = 'MS_batch',
                       facet_col = 'instrument',
                       color_by_batch = T, color_scheme = 'brewer',
                       theme = 'classic',
                       plot_title = NULL, order_per_facet = F){
  
  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(df_long[[sample_id_col]])) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  
  if (!all(c(batch_col, sample_id_col) %in% names(df_long))){
    if (!is.null(sample_annotation)){
      df_long = df_long %>% merge(sample_annotation,
                                            by = sample_id_col)
      if(is.numeric(df_long[[batch_col]])){
        df_long[,batch_col] <- as.factor(df_long[, batch_col])
      }
    } else {
      if (color_by_batch){
        stop('batches cannot be colored if the batch column cannot be defined,
             check sample_annotation and data matrix')
      }
    }
  }
  if (is.null(order_col)){
    warning('order column not defined, taking order of files in the data matrix instead')
    order_col = 'order_temp_col'
    df_long[[order_col]] = match(df_long[[sample_id_col]],
                                         unique(df_long[[sample_id_col]]))
  } else if (!(order_col %in% names(sample_annotation)) &
             !(order_col %in% names(df_long))){
    warning('order column not found in sample annotation, taking order of files in the data matrix instead')
    order_col = 'order_temp_col'
    df_long[[order_col]] = match(df_long[[sample_id_col]],
                                         unique(df_long[[sample_id_col]]))
    order_per_facet = T
  }

  if (order_per_facet){
    if (!is.null(facet_col)){
      warning('defining order within each facet')
      df_long = df_long %>%
        group_by_at(vars(one_of(facet_col))) %>%
        mutate(order = rank(UQ(sym(order_col))))
    }
  }

  gg = ggplot(df_long, aes_string(x = order_col, y = measure_col,
                                       group = order_col))+
    geom_boxplot(outlier.size = .15)+
    theme_bw()+
    theme(plot.title = element_text(face = 'bold', hjust = .5))


  if (!is.numeric(df_long[[order_col]])){
    warning(sprintf('order column is not numeric, assuming the order is irrelevant
                    of %s order follows the run order.', order_col))
    df_long[[order_col]] = factor(df_long[[order_col]],
                                              levels = df_long[[order_col]])
    gg = gg +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  }
  if(color_by_batch){
    gg = gg + aes_string(fill = batch_col)
    if(length(color_scheme) == 1 & color_scheme == 'brewer'){
      n_batches <- length(unique(df_long[[batch_col]]))
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
  if(!is.null(facet_col)){
    if(!(facet_col %in% names(df_long))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present in the data,
                   check sample annotation data frame', facet_col))
    }
    gg = gg + facet_wrap(as.formula(paste("~", facet_col)),
                         dir = 'v', scales = "free_x")
  }
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)+
      theme(plot.title = element_text(hjust = .5, face = 'bold', size = 16))
  }
  if (max(df_long[[order_col]]) > 100){
    gg = gg + theme(legend.position="top")
  }
  return(gg)
}


#' cluster the data matrix to visually inspect which confounder dominates
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param color_df data frame of colors, as created by
#'   `sample_annotation_to_colors`
#' @param distance distance metric used for clustering
#' @param agglomeration agglomeration methods as used by `hclust`
#' @param label_samples if \code{TRUE} sample IDs (column names of \code{data_matrix}) will be printed
#' @param label_font size of the font. Is active if \code{label_samples} is \code{TRUE}, ignored otherwise
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param ... other parameters of `plotDendroAndColors` from `WGCNA` package
#'
#' @export
#'
#' @examples
#' @seealso \code{\link[stats]{hclust}}, \code{\link{sample_annotation_to_colors}},
#'   \code{\link[WGCNA]{plotDendroAndColors}}
plot_hierarchical_clustering  <- function(data_matrix, color_df,
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

#' Plot the heatmap of samples
#'
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. in most
#'   function, it is assumed that this is the log transformed version of the
#'   original data
#' @param sample_annotation data matrix with \enumerate{
#' \item  `sample_id_col` (this can be repeated as row names)
#'   \item  biological and
#'   \item  technical covariates (batches etc)
#' }; each column of sample annotation will get it's own row. If `cluster_cols = T` this will indicate,
#' whether sample proximity is driven by one of biolical or technical factors
#' @param fill_the_missing boolean value determining if missing values should be
#'   substituted with -1 (and colored with black)
#' @param cluster_rows boolean value determining if rows should be clustered
#' @param cluster_cols boolean value determining if columns should be clustered
#' @param sample_annotation_col biological or technical factors to be annotated in heatmap columns
#' @param sample_annotation_row biological or technical factors to be annotated in heatmap rows 
#' @param annotation_color_list list specifying colors for columns (samples).
#'   Best created by `sample_annotation_to_colors`
#' @param heatmap_color vector of colors used in heatmap (typicall a gradient)
#' @param color_for_missing special color to make missing values. Usually black or white, depending on `heatmap_color`
#' @param filename filepath where to save the image
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param ... other parameters of \code{link[pheatmap]{pheatmap}}
#' 
#' @return object returned by \code{link[pheatmap]{pheatmap}}
#' @export
#'
#' @examples
#' @seealso \code{\link{sample_annotation_to_colors}}, \code{\link[pheatmap]{pheatmap}}
plot_heatmap <- function(data_matrix, sample_annotation = NULL, sample_id_col = 'FullRunName',
                         sample_annotation_col = NULL, 
                         sample_annotation_row = NULL, 
                         fill_the_missing = T, cluster_rows = T, cluster_cols = F,
                         annotation_color_list = NA,
                         heatmap_color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
                         color_for_missing = 'black',
                         filename = NA, plot_title = NA,
                         ...){
  if(fill_the_missing) {
    data_matrix[is.na(data_matrix)] = 0
    warning('substituting missing values with 0, this might affect the clustering!')
    heatmap_color = c(color_for_missing, heatmap_color)
  }
  
  if (is.null(sample_annotation)){
    annotation_col = NA
    annotation_row = NA
  }
  
  if(!is.null(sample_annotation)){
    if(!is.null(sample_annotation_col) && is.null(sample_annotation_row)){
      annotation_row = NA
      annotation_col = sample_annotation %>% 
        select(one_of(sample_id_col, sample_annotation_col)) %>%
        remove_rownames %>% 
        column_to_rownames(var=sample_id_col)  
      
    }
    if(!is.null(sample_annotation_row) && is.null(sample_annotatation_col)){
      annotation_col = NA
      annotation_row = sample_annotation %>% 
        select(one_of(sample_id_col, sample_annotation_row)) %>%
        remove_rownames %>% 
        column_to_rownames(var=sample_id_col)
    }
    if(is.null(sample_annotation_col) && is.null(sample_annotation_row)){
      warning("Sample annotation columns and rows are not specified for heatmap.")
      annotation_col = sample_annotation %>% 
        remove_rownames %>% 
        column_to_rownames(var=sample_id_col)
    }
  }
  
  p <- pheatmap(data_matrix, cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                color = heatmap_color,
                annotation_col = annotation_col, annotation_row = annotation_row, 
                annotation_colors = annotation_color_list,
                filename = filename, main = plot_title, ...)
  return(p)
}


calculate_PVCA <- function(data_matrix, sample_annotation, factors_for_PVCA,
                 threshold_pca, threshold_var = Inf) {

  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(colnames(data_matrix))) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  
  covrts.annodf = Biobase::AnnotatedDataFrame(data=sample_annotation)
  expr_set = Biobase::ExpressionSet(data_matrix[,rownames(sample_annotation)], covrts.annodf)
  pvcaAssess = pvca::pvcaBatchAssess (expr_set, factors_for_PVCA, threshold = threshold_pca)
  pvcaAssess_df = data.frame(weights = as.vector(pvcaAssess$dat),
                             label = pvcaAssess$label,
                             stringsAsFactors = F)

  label_of_small = sprintf('Below %1.0f%%', 100*threshold_var)
  if (sum(pvcaAssess_df$weights < threshold_var) > 1){
    pvca_res_small = sum(pvcaAssess_df$weights[pvcaAssess_df$weights < threshold_var])
    pvca_res = pvcaAssess_df[pvcaAssess_df$weights >= threshold_var, ]
    pvca_res_add = data.frame(weights = pvca_res_small, label = label_of_small)
    pvca_res = rbind(pvca_res, pvca_res_add)
  } else {
    pvca_res = pvcaAssess_df
    }
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
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param feature_id_col name of the column with feature/gene/peptide/protein
#'   ID used in the long format representation \code{df_long}. In the wide
#'   formatted representation \code{data_matrix} this corresponds to the row
#'   names.
#' @param technical_covariates vector `sample_annotation` column names that are
#'   technical covariates
#' @param biological_covariates vector `sample_annotation` column names, that
#'   are biologically meaningful covariates
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#' @param threshold_pca the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param threshold_var the percentile value of weight each of the covariates
#'   needs to explain (the rest will be lumped together)
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param feature_id_col name of the column with feature/gene/peptide/protein
#'   ID used in the long format representation \code{df_long}. In the wide
#'   formatted representation \code{data_matrix} this corresponds to the row
#'   names.
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see
#'   examples)
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#'
#' @return list of two items: plot =gg, df = pvca_res
#' @export
#'
#' @examples
#' @seealso \code{\link{sample_annotation_to_colors}}, \code{\link[ggplot2]{ggplot}}
plot_PVCA <- function(data_matrix, sample_annotation,
                      sample_id_col = 'FullRunName',
                      feature_id_col = 'peptide_group_label',
                      technical_covariates = c('MS_batch', 'instrument'),
                      biological_covariates = c('cell_line','drug_dose'),
                      fill_the_missing = 0,
                      threshold_pca = .6, threshold_var = .01,
                      colors_for_bars = NULL,
                      theme = 'classic', plot_title = NULL){
  factors_for_PVCA = c(technical_covariates, biological_covariates)

  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(colnames(data_matrix))) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  
  sample_names = sample_annotation[[sample_id_col]]
  sample_annotation = sample_annotation %>% select(one_of(factors_for_PVCA)) %>%
    mutate_if(lubridate::is.POSIXct, as.numeric)
  sample_annotation = as.data.frame(sample_annotation)
  rownames(sample_annotation) = sample_names

  if (!is.null(sample_id_col)){
    if(sample_id_col %in% names(sample_annotation)){
      rownames(sample_annotation) = sample_annotation[[sample_id_col]]
    }
    if(feature_id_col %in% colnames(data_matrix)){
      if(is.data.frame(data_matrix)){
        warning(sprintf('feature_id_col with name %s in data matrix instead of rownames,
                        this might cause errors in other diagnostic functions,
                        assign values of this column to rowname and remove from the data frame!', feature_id_col))
        rownames(data_matrix) = data_matrix[[feature_id_col]]
        data_matrix[[feature_id_col]] = NULL
        data_matrix = as.matrix(data_matrix)
        }

      }
  } else {
      if (!all(rownames(sample_annotation) %in% colnames(data_matrix))){
        stop('sample names differ between data matrix and sample annotation')
      }
  }

  if (any(is.na(as.vector(data_matrix)))){
    warning('PVCA cannot operate with missing values in the matrix')
    if(!is.null(fill_the_missing)){
      warning(sprintf('filling missing value with %s', fill_the_missing))
      data_matrix[is.na(data_matrix)] = fill_the_missing
    } else {
      warning('filling value is NULL, removing features with missing values')
      data_matrix = data_matrix[complete.cases(data_matrix),]
    }
  }

  pvca_res = calculate_PVCA(data_matrix, sample_annotation, factors_for_PVCA,
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

  label_of_small = sprintf('Below %1.0f%%', 100*threshold_var)
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
  pvca_res = pvca_res %>%
    mutate(label = factor(label, levels=label))

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

  } else {
    if (length(colors_for_bars) != 4){
      color_names = paste(c('residual', 'biological', 'biol:techn', 'technical'), collapse = ' ')
      warning(sprintf('four colors for: %s were expected', color_names))
    }
  }
  gg = gg + scale_fill_manual(values = colors_for_bars)

  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
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
#' @param feature_id_col name of the column with feature/gene/peptide/protein
#'   ID used in the long format representation \code{df_long}. In the wide
#'   formatted representation \code{data_matrix} this corresponds to the row
#'   names.
#' @param color_by column name (as in `sample_annotation`) to color by
#' @param PC_to_plot principal component numbers for x and y axis
#' @param colors_for_factor named vector of colors for the `color_by` variable
#' @param fill_the_missing boolean value determining if missing values should be
#'   substituted with -1 (and colored with black)
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see
#'   examples)
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#'
#' @return ggplot scatterplot colored by factor levels of column specified in
#'   `factor_to_color`
#' @export
#'
#' @examples
#' @seealso \code{\link[ggfortify]{autoplot.pca_common}}, \code{\link[ggplot2]{ggplot}}
plot_PCA <- function(data_matrix, sample_annotation,
                     feature_id_col = 'peptide_group_label',
                     color_by = 'MS_batch',
                     PC_to_plot = c(1,2), fill_the_missing = 0,
                     colors_for_factor = NULL,
                     theme = 'classic',
                     plot_title = NULL){

  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(colnames(data_matrix))) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  
  if(!is.null(feature_id_col)){
    if(feature_id_col %in% colnames(data_matrix)){
      if(is.data.frame(data_matrix)){
        warning(sprintf('feature_id_col with name %s in data matrix instead of rownames,
                        this might cause errors in other diagnostic functions,
                        assign values of this column to rowname and remove from the data frame!', feature_id_col))
      }
      rownames(data_matrix) = data_matrix[[feature_id_col]]
      data_matrix[[feature_id_col]] = NULL
      data_matrix = as.matrix(data_matrix)
      }
  }

  if (any(is.na(as.vector(data_matrix)))){
    warning('PCA cannot operate with missing values in the matrix')
    if(!is.null(fill_the_missing)){
      if (!is.numeric(fill_the_missing)){
        fill_the_missing = 0
      }
      warning(sprintf('filling missing value with %s', fill_the_missing))
      data_matrix[is.na(data_matrix)] = fill_the_missing
    } else {
      warning('filling value is NULL, removing features with missing values')
      data_matrix = data_matrix[complete.cases(data_matrix),]
    }
  }

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
    gg = gg + aes_string(color = color_by) + scale_color_manual(values = colors_for_factor)
  }
  if(!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  return(gg)
}


