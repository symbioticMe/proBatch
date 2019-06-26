#' cluster the data matrix to visually inspect which confounder dominates
#'
#' @inheritParams proBatch
#' @param color_df data frame of colors, as created by 
#' \code{sample_annotation_to_colors}
#' @param distance distance metric used for clustering
#' @param agglomeration agglomeration methods as used by \code{hclust}
#' @param label_samples if \code{TRUE} sample IDs (column names of 
#' \code{data_matrix}) will be printed
#' @param label_font size of the font. Is active if \code{label_samples} is 
#' \code{TRUE}, ignored otherwise
#' @param fill_the_missing boolean value determining if  missing values 
#' should be substituted with -1 (and colored with   \code{color_for_missing}). 
#' If \code{NULL}, features with missing values are excluded.
#' @param ... other parameters of \code{plotDendroAndColors} from \code{WGCNA} package
#'
#' @return No return
#' @examples
#' color_scheme <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch','EarTag', "Strain", "Diet", "digestion_batch", "Sex"),
#' date_columns = 'DateTime',
#' numeric_columns = c('order'))
#' 
#' color_annotation <- color_scheme$color_df
#' 
#' hiarchical_clustering_plot <- plot_hierarchical_clustering(
#' example_proteome_matrix, color_annotation,  
#' distance = "euclidean", agglomeration = 'complete',
#' label_samples = FALSE)
#' 
#' @export
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link{sample_annotation_to_colors}},
#'   \code{\link[WGCNA]{plotDendroAndColors}}
plot_hierarchical_clustering  <- function(data_matrix, color_df,
                                          fill_the_missing = 0,
                                          distance = "euclidean",
                                          agglomeration = 'complete',
                                          label_samples = TRUE, label_font = .2,
                                          filename = NULL, width = 38, height = 25, 
                                          units = c('cm','in','mm'), 
                                          plot_title = NULL,
                                          ...){
  if (any(is.na(as.vector(data_matrix)))){
    warning('Hierarchical clustering cannot operate with missing values in the matrix')
    if(!is.null(fill_the_missing)){
      if (!is.numeric(fill_the_missing)){
        fill_the_missing = 0
      } else{
        warning(sprintf('filling missing value with %s', fill_the_missing))
        data_matrix[is.na(data_matrix)] = fill_the_missing
      }
    } else {
      warning('filling value is NULL, removing features with missing values')
      data_matrix = data_matrix[complete.cases(data_matrix),]
    }
  }
  
  dist_matrix = dist(t(as.matrix(data_matrix)), method = distance)
  hierarchical_clust = hclust(dist_matrix, method = agglomeration)
  if (label_samples){
    cex.dendroLabels = label_font
    if (ncol(data_matrix) > 80){
      warning('Too many samples, adjust the font with `label_font` argument or
              remove labels by setting `label_samples = FALSE` in function call')
    }
  } else{
    cex.dendroLabels = 0.9
  }
  
  if (is.null(filename)){
    plotDendroAndColors(hierarchical_clust, color_df, rowTextAlignment = 'left',
                        main = plot_title,
                        hang = -0.1, addGuide = TRUE, dendroLabels = FALSE, 
                        cex.dendroLabels = cex.dendroLabels,
                        ...)
  } else {
    
    units_adjusted = adjust_units(units, width, height)
    units = units_adjusted$unit
    width = units_adjusted$width
    height = units_adjusted$height
    
    if (is.na(width)){
      width = 7
    }
    if (is.na(height)){
      height = 7
    }
    if (file_ext(filename) == 'pdf'){
      pdf(file = filename, width = width, height = height, title = plot_title)
    } else if(file_ext(filename) == 'png'){
      png(file = filename, width = width, height = height, units = units, res = 300)
    } else{
      stop('currently only pdf and png extensions for filename are implemented')
    }
    
    plotDendroAndColors(hierarchical_clust, color_df, rowTextAlignment = 'left',
                        main = plot_title,
                        hang = -0.1, addGuide = TRUE, dendroLabels = FALSE, 
                        cex.dendroLabels = cex.dendroLabels,
                        ...)
    
    dev.off()
  }
}

#' Plot the heatmap of samples
#'
#' @inheritParams proBatch
#' @param sample_annotation data matrix with \enumerate{
#' \item  \code{sample_id_col} (this can be repeated as row names)
#'   \item  biological and
#'   \item  technical covariates (batches etc)
#' }; each column of sample annotation will get it's own row. 
#' If \code{cluster_cols = T} this will indicate,
#' whether sample proximity is driven by one of 
#' biological or technical factors
#' @param sample_id_col name of the column in 
#' sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param fill_the_missing numeric value that the missing values are
#'   substituted with, or \code{NULL} if features with missing values are to be excluded.
#' @param cluster_rows boolean value determining if rows 
#' should be clustered
#' @param cluster_cols boolean value determining if columns 
#' should be clustered
#' @param sample_annotation_col biological or technical 
#' factors to be annotated in heatmap columns
#' @param sample_annotation_row biological or technical 
#' factors to be annotated in heatmap rows 
#' @param annotation_color_list list specifying colors 
#' for columns (samples). Best created by \code{sample_annotation_to_colors}
#' @param heatmap_color vector of colors used in heatmap (typicall a gradient)
#' @param color_for_missing special color to make missing values. 
#' Usually black or white, depending on \code{heatmap_color}
#' @param ... other parameters of \code{link[pheatmap]{pheatmap}}
#' 
#' @return object returned by \code{link[pheatmap]{pheatmap}}
#' @export
#' 
#' @examples 
#' color_scheme <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch','EarTag', "Strain", 
#' "Diet", "digestion_batch", "Sex"),
#' date_columns = 'DateTime',
#' numeric_columns = c('order'))
#' 
#' heatmap_plot <- plot_heatmap(log_transform_dm(example_proteome_matrix), 
#' example_sample_annotation, 
#' sample_annotation_col = c("MS_batch",  "digestion_batch", "Diet"), 
#' cluster_cols = TRUE, 
#' annotation_color_list = color_scheme$list_of_colors,
#' show_rownames = FALSE, show_colnames = FALSE)
#'
#' @seealso \code{\link{sample_annotation_to_colors}}, \code{\link[pheatmap]{pheatmap}}
plot_heatmap <- function(data_matrix, sample_annotation = NULL, sample_id_col = 'FullRunName',
                         sample_annotation_col = NULL, 
                         sample_annotation_row = NULL, 
                         fill_the_missing = -1, cluster_rows = TRUE, cluster_cols = FALSE,
                         annotation_color_list = NA,
                         heatmap_color = colorRampPalette(
                           rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                         color_for_missing = 'black',
                         filename = NULL, width = 7, height = 7, 
                         units = c('cm','in','mm'), 
                         plot_title = NULL,
                         ...){
  
  df_long = matrix_to_long(data_matrix, sample_id_col = sample_id_col)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     merge = F)
  data_matrix = long_to_matrix(df_long, sample_id_col = sample_id_col)
  rm(df_long)
  
  if (any(is.na(as.vector(data_matrix)))){
    warning('Heatmap cannot operate with missing values in the matrix')
    if(!is.null(fill_the_missing)){
      heatmap_color = c(color_for_missing, heatmap_color)
      if (!is.numeric(fill_the_missing)){
        fill_the_missing = 0
      } else{
        warning(sprintf('filling missing value with %s', fill_the_missing))
        data_matrix[is.na(data_matrix)] = fill_the_missing
      }
    } else {
      warning('filling value is NULL, removing features with missing values')
      data_matrix = data_matrix[complete.cases(data_matrix),]
    }
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
    if(!is.null(sample_annotation_row) && is.null(sample_annotation_col)){
      annotation_col = NA
      annotation_row = sample_annotation %>% 
        select(one_of(sample_id_col, sample_annotation_row)) %>%
        remove_rownames %>% 
        column_to_rownames(var=sample_id_col)
    }
    if(is.null(sample_annotation_col) && is.null(sample_annotation_row)){
      warning("annotation_row and annotation_col are not specified for heatmap, 
              using sample_annotation for annotation_col.")
      annotation_col = sample_annotation %>% 
        remove_rownames %>% 
        column_to_rownames(var=sample_id_col)
      annotation_row = NA
    }
  }
  
  if(is.null(plot_title)){
    plot_title = NA
  }
  if(is.null(filename)){
    filename = NA
  }
  units_adjusted = adjust_units(units, width, height)
  units = units_adjusted$unit
  width = units_adjusted$width
  height = units_adjusted$height
  p <- pheatmap(data_matrix, 
                cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                color = heatmap_color,
                annotation_col = annotation_col, annotation_row = annotation_row, 
                annotation_colors = annotation_color_list,
                filename = filename, width = width, height = height,
                main = plot_title, ...)
  return(p)
}


calculate_PVCA <- function(data_matrix, sample_annotation, factors_for_PVCA,
                           pca_threshold, variance_threshold = Inf) {
  
  covrts.annodf = Biobase::AnnotatedDataFrame(data=sample_annotation)
  expr_set = Biobase::ExpressionSet(data_matrix, covrts.annodf)
  pvcaAssess = pvcaBatchAssess(expr_set, factors_for_PVCA, threshold = pca_threshold)
  pvcaAssess_df = data.frame(weights = as.vector(pvcaAssess$dat),
                             label = pvcaAssess$label,
                             stringsAsFactors = FALSE)
  
  label_of_small = sprintf('Below %1.0f%%', 100*variance_threshold)
  if (sum(pvcaAssess_df$weights < variance_threshold) > 1){
    pvca_res_small = sum(pvcaAssess_df$weights[pvcaAssess_df$weights < variance_threshold])
    pvca_res = pvcaAssess_df[pvcaAssess_df$weights >= variance_threshold, ]
    pvca_res_add = data.frame(weights = pvca_res_small, label = label_of_small)
    pvca_res = rbind(pvca_res, pvca_res_add)
  } else {
    pvca_res = pvcaAssess_df
  }
  return(pvca_res)
}

#' Plot variance distribution by variable
#' 
#' @inheritParams proBatch
#' @param technical_covariates vector \code{sample_annotation} column names that are
#'   technical covariates
#' @param biological_covariates vector \code{sample_annotation} column names, that
#'   are biologically meaningful covariates
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#' @param pca_threshold the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param variance_threshold the percentile value of weight each of the covariates
#'   needs to explain (the rest will be lumped together)
#' @param fill_the_missing boolean value determining if  missing values 
#' should be substituted with -1 (and colored with   \code{color_for_missing}). 
#' If \code{NULL}, features with missing values are excluded.
#'
#' @return list of two items: plot =gg, df = pvca_res
#' @export
#'
#' @examples 
#' matrix_test <- example_proteome_matrix[1:50, ]
#' pvca_plot <- plot_PVCA(matrix_test, example_sample_annotation, 
#' technical_covariates = c('MS_batch', 'digestion_batch'),
#' biological_covariates = c("Diet", "Sex", "Strain"))
#' 
#' \dontrun{
#' pvca_plot <- plot_PVCA(matrix_test, example_sample_annotation, 
#' technical_covariates = c('MS_batch', 'digestion_batch'),
#' biological_covariates = c("Diet", "Sex", "Strain"), 
#' filename = 'test_PVCA.png', width = 28, height = 22, units = 'cm')
#' }
#' 
#' @seealso \code{\link{sample_annotation_to_colors}}, 
#' \code{\link[ggplot2]{ggplot}}
plot_PVCA <- function(data_matrix, sample_annotation,
                      sample_id_col = 'FullRunName',
                      feature_id_col = 'peptide_group_label',
                      technical_covariates = c('MS_batch', 'instrument'),
                      biological_covariates = c('cell_line','drug_dose'),
                      fill_the_missing = -1,
                      pca_threshold = .6, variance_threshold = .01,
                      colors_for_bars = NULL,
                      filename = NULL, width = NA, height = NA, 
                      units = c('cm','in','mm'),
                      plot_title = NULL,
                      theme = 'classic'){
  factors_for_PVCA = c(technical_covariates, biological_covariates)
  
  df_long = matrix_to_long(data_matrix, sample_id_col = sample_id_col)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col = NULL, 
                                     facet_col = NULL, merge = FALSE)
  data_matrix = long_to_matrix(df_long, sample_id_col = sample_id_col)
  
  sample_annotation = sample_annotation   %>% 
    select(one_of(c(sample_id_col, factors_for_PVCA))) %>%
    mutate_if(is.POSIXct, as.numeric) %>%
    as.data.frame()%>%
    column_to_rownames(var = sample_id_col)
  
  data_matrix = check_feature_id_col_in_dm(feature_id_col, data_matrix)
  
  if (any(is.na(as.vector(data_matrix)))){
    warning('PVCA cannot operate with missing values in the matrix')
    if(!is.null(fill_the_missing)){
      if (!is.numeric(fill_the_missing)){
        fill_the_missing = 0
      } else{
        warning(sprintf('filling missing value with %s', fill_the_missing))
        data_matrix[is.na(data_matrix)] = fill_the_missing
      }
    } else {
      warning('filling value is NULL, removing features with missing values')
      data_matrix = data_matrix[complete.cases(data_matrix),]
    }
  }
  
  pvca_res = calculate_PVCA(data_matrix, sample_annotation, factors_for_PVCA,
                            pca_threshold, variance_threshold = variance_threshold)
  
  tech_interactions = expand.grid(technical_covariates, 
                                  technical_covariates) %>%
    mutate(tech_interactions = paste(Var1, Var2, sep = ':')) %>%
    pull(tech_interactions)
  biol_interactions = expand.grid(biological_covariates, 
                                  biological_covariates) %>%
    mutate(biol_interactions = paste(Var1, Var2, sep = ':')) %>%
    pull(biol_interactions)
  
  if(!(all(rownames(sample_annotation) %in% colnames(data_matrix)))){
    stop('check data matrix column names or these in sample annotation')
  }
  
  label_of_small = sprintf('Below %1.0f%%', 100*variance_threshold)
  technical_covariates = c(technical_covariates, tech_interactions)
  biological_covariates = c(biological_covariates, biol_interactions)
  pvca_res = pvca_res %>% mutate(category = ifelse(label %in% technical_covariates, 
                                                   'technical',
                                                   ifelse(label %in% biological_covariates, 
                                                          'biological',
                                                          ifelse(label %in% c(label_of_small, 'resid'), 
                                                                 'residual', 'biol:techn'))))
  
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
  
  
  if(is.null(colors_for_bars)){
    colors_for_bars = c('grey', wes_palettes$Rushmore[3:5])
    names(colors_for_bars) = c('residual', 'biological', 
                               'biol:techn', 'technical')
    
  } else {
    if (length(colors_for_bars) != 4){
      color_names = paste(c('residual', 'biological', 'biol:techn', 
                            'technical'), collapse = ' ')
      warning(sprintf('four colors for: %s were expected', color_names))
    }
  }
  gg = gg + scale_fill_manual(values = colors_for_bars)
  
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  
  #Change the theme
  if(!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic()
  }else{
    message("plotting with default ggplot theme, only theme = 'classic' implemented")
  }
  
  save_ggplot(filename, units, width, height, gg)
  
  gg = gg +
    theme(axis.title.x = NULL, 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
    xlab(NULL)+
    theme(text = element_text(size=15))+
    guides(fill=guide_legend(override.aes=list(color=NA), title=NULL))
  
  return(list(plot =gg, df = pvca_res))
}

#' plot PCA plot
#'
#' @inheritParams proBatch
#' @param color_by column name (as in \code{sample_annotation}) to color by
#' @param PC_to_plot principal component numbers for x and y axis
#' @param fill_the_missing boolean value determining if  missing values 
#' should be substituted with -1 (and colored with   \code{color_for_missing}). 
#' If \code{NULL}, features with missing values are excluded.
#'
#' @return ggplot scatterplot colored by factor levels of column specified in
#'   \code{factor_to_color}
#' @export
#'
#' @examples 
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation, 
#' color_by = 'MS_batch', plot_title = "PCA colored by MS batch")
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation, 
#' color_by = 'digestion_batch', plot_title = "PCA colored by digestion_batch")
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation, 
#' color_by = 'DateTime', plot_title = "PCA colored by DateTime")
#' 
#' \dontrun{
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation, 
#' color_by = 'DateTime', plot_title = "PCA colored by DateTime",
#' filename = 'test_PCA.png', width = 14, height = 9, units = 'cm')
#' }
#' 
#' @seealso \code{\link[ggfortify]{autoplot.pca_common}}, 
#' \code{\link[ggplot2]{ggplot}}
plot_PCA <- function(data_matrix, sample_annotation,
                     feature_id_col = 'peptide_group_label',
                     color_by = 'MS_batch',
                     PC_to_plot = c(1,2), fill_the_missing = -1,
                     color_scheme = 'brewer',
                     filename = NULL, width = NA, height = NA, 
                     units = c('cm','in','mm'),
                     plot_title = NULL,
                     theme = 'classic'){
  
  data_matrix = check_feature_id_col_in_dm(feature_id_col, data_matrix)
  
  if (any(is.na(as.vector(data_matrix)))){
    warning('PCA cannot operate with missing values in the matrix')
    if(!is.null(fill_the_missing)){
      if (!is.numeric(fill_the_missing)){
        fill_the_missing = 0
      } else{
        warning(sprintf('filling missing value with %s', fill_the_missing))
        data_matrix[is.na(data_matrix)] = fill_the_missing
      }
    } else {
      warning('filling value is NULL, removing features with missing values')
      data_matrix = data_matrix[complete.cases(data_matrix),]
    }
  }
  
  
  pr_comp_res <- prcomp(t(data_matrix))
  gg = autoplot(pr_comp_res, data = sample_annotation,
                x = PC_to_plot[1], y = PC_to_plot[2])
  
  
  #add colors
  
  if(length(color_by) > 1){
    warning('Coloring by the first column specified')
    color_by = color_by[1]
  } #TODO: create the ggpubr graph with multiple panels, colored by factors
  gg = color_by_factor(color_by_batch = TRUE, 
                       batch_col = color_by, gg = gg, 
                       color_scheme = color_scheme, 
                       sample_annotation = sample_annotation,
                       fill_or_color = 'color')
  
  #Add the title
  if(!is.null(plot_title)) {
    gg = gg + ggtitle(plot_title)+
      theme(plot.title = element_text(face = 'bold',hjust = .5))
  }
  
  #Change the theme
  if(!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic()
  }else{
    message("plotting with default ggplot theme, only theme = 'classic' implemented")
  }
  
  save_ggplot(filename, units, width, height, gg)
  
  return(gg)
}


