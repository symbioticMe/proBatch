#' cluster the data matrix to visually inspect which confounder dominates
#'
#' @inheritParams proBatch
#' @param distance distance metric used for clustering
#' @param agglomeration agglomeration methods as used by \code{hclust}
#' @param label_samples if \code{TRUE} sample IDs (column names of 
#' \code{data_matrix}) will be printed
#' @param label_font size of the font. Is active if \code{label_samples} is 
#' \code{TRUE}, ignored otherwise
#' @param fill_the_missing numeric value determining how  missing values 
#' should be substituted. If \code{NULL}, features with missing values are 
#' excluded.
#' @param ... other parameters of \code{plotDendroAndColors} from 
#' \code{WGCNA} package
#'
#' @return No return
#' @examples
#' 
#' selected_batches = example_sample_annotation$MS_batch %in% 
#'                                               c('Batch_1', 'Batch_2')
#' selected_samples = example_sample_annotation$FullRunName[selected_batches]
#' test_matrix = example_proteome_matrix[,selected_samples]
#' 
#' hierarchical_clustering_plot <- plot_hierarchical_clustering(
#' example_proteome_matrix, example_sample_annotation,
#' factors_to_plot = c('MS_batch', 'Diet', 'DateTime'),
#' color_list = NULL,  
#' distance = "euclidean", agglomeration = 'complete',
#' label_samples = FALSE)
#' 
#' #with defined color scheme:
#' color_list <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch', "Strain", "Diet", "digestion_batch"),
#' numeric_columns = c('DateTime', 'order'))
#' hierarchical_clustering_plot <- plot_hierarchical_clustering(
#' example_proteome_matrix, example_sample_annotation,
#' factors_to_plot = c('MS_batch', "Strain", 'DateTime', "digestion_batch"),
#' color_list = color_list,  
#' distance = "euclidean", agglomeration = 'complete',
#' label_samples = FALSE)
#' 
#' @export
#'
#' @seealso \code{\link[stats]{hclust}}, 
#'   \code{\link{sample_annotation_to_colors}},
#'   \code{\link[WGCNA]{plotDendroAndColors}}
plot_hierarchical_clustering  <- function(data_matrix, sample_annotation,
                                          sample_id_col = 'FullRunName',
                                          color_list = NULL,
                                          factors_to_plot = NULL,
                                          fill_the_missing = 0,
                                          distance = "euclidean",
                                          agglomeration = 'complete',
                                          label_samples = TRUE, label_font = .2,
                                          filename = NULL, 
                                          width = 38, height = 25, 
                                          units = c('cm','in','mm'), 
                                          plot_title = NULL,
                                          ...){

  df_long = matrix_to_long(data_matrix, sample_id_col = sample_id_col)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     merge = FALSE)
  data_matrix = long_to_matrix(df_long, sample_id_col = sample_id_col)
  rm(df_long)
  
  warning_message <- 'Hierarchical clustering cannot operate with missing values
                      in the matrix'
  data_matrix = handle_missing_values(data_matrix, warning_message, 
                                      fill_the_missing)
  
  dist_matrix = dist(t(as.matrix(data_matrix)), method = distance)
  hierarchical_clust = hclust(dist_matrix, method = agglomeration)
  if (label_samples){
    cex.dendroLabels = label_font
    if (ncol(data_matrix) > 80){
      warning('Too many samples, adjust the font with `label_font` argument or
              remove labels by setting `label_samples = FALSE` in 
              function call')
    }
  } else{
    cex.dendroLabels = 0.9
  }
  
  factors_without_colors = setdiff(factors_to_plot, names(color_list))
  if(length(factors_without_colors) > 0){
    warning('color_list for samples annotation not defined, inferring 
             automatically. Numeric/factor columns are guessed, for more 
            controlled color mapping use sample_annotation_to_colors()')
    color_list_new <- sample_annotation_to_colors(sample_annotation, 
                                                sample_id_col = sample_id_col,
                                                factor_columns = factors_without_colors, 
                                                numeric_columns = NULL)
    color_list = c(color_list, color_list_new)
  }
  
  
  if (length(setdiff(names(color_list), factors_to_plot)) > 0 && 
      !is.null(factors_to_plot)){
    color_list = color_list[factors_to_plot]
  }
  
  #transform color list to color_df
  sample_annotation = sample_annotation %>% 
    filter(!!sym(sample_id_col) %in% colnames(data_matrix))
  color_df = color_list_to_df(color_list, sample_annotation, sample_id_col)
  
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
      png(filename = filename, width = width, height = height, units = units, 
          res = 300)
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

#' Plot the heatmap of samples (cols) vs features (rows)
#'
#' @inheritParams proBatch
#' @param factors_to_plot vector of technical and biological factors to be 
#' plotted in this diagnostic plot (assumed to be present in 
#' \code{sample_annotation})
#' @param fill_the_missing numeric value that the missing values are
#'   substituted with, or \code{NULL} if features with missing values are to be 
#'   excluded.
#' @param color_for_missing special color to make missing values. 
#' Usually black or white, depending on \code{heatmap_color}
#' @param heatmap_color vector of colors used in heatmap (typicall a gradient)
#' @param cluster_rows boolean value determining if rows 
#' should be clustered
#' @param cluster_cols boolean value determining if columns 
#' should be clustered
#' @param factors_of_feature_ann vector of factors that characterize features, 
#' as listed in \code{peptide_annotation}
#' @param color_list_features list, as returned by 
#' \code{sample_annotation_to_colors}, 
#' but mapping \code{peptide_annotation} where each item contains a color vector 
#' for each factor to be mapped to the color.
#' @param ... other parameters of \code{link[pheatmap]{pheatmap}}
#' 
#' @return object returned by \code{link[pheatmap]{pheatmap}}
#' @export
#' 
#' @examples 
#' 
#' log_transformed_matrix = log_transform_dm(example_proteome_matrix)
#' heatmap_plot <- plot_heatmap_diagnostic(log_transformed_matrix, 
#' example_sample_annotation, 
#' factors_to_plot = c("MS_batch",  "digestion_batch", "Diet", 'DateTime'), 
#' cluster_cols = TRUE, cluster_rows = FALSE,
#' show_rownames = FALSE, show_colnames = FALSE)
#' 
#' color_list <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch','EarTag', "Strain", 
#' "Diet", "digestion_batch", "Sex"),
#' numeric_columns = c('DateTime', 'order'))
#' 
#' log_transformed_matrix = log_transform_dm(example_proteome_matrix)
#' heatmap_plot <- plot_heatmap_diagnostic(log_transformed_matrix, 
#' example_sample_annotation, 
#' factors_to_plot = c("MS_batch",  "digestion_batch", "Diet", 'DateTime'), 
#' cluster_cols = TRUE, cluster_rows = FALSE,
#' color_list = color_list,
#' show_rownames = FALSE, show_colnames = FALSE)
#'
#' @seealso \code{\link{sample_annotation_to_colors}}, 
#' \code{\link[pheatmap]{pheatmap}}
#' 
#' @name plot_heatmap_diagnostic
plot_heatmap_diagnostic <- function(data_matrix, sample_annotation = NULL, 
                                    sample_id_col = 'FullRunName',
                                    factors_to_plot = NULL,
                                    fill_the_missing = -1, 
                                    color_for_missing = 'black',
                                    heatmap_color = colorRampPalette(
                                      rev(brewer.pal(n = 7, 
                                                     name = "RdYlBu")))(100),
                                    cluster_rows = TRUE, cluster_cols = FALSE,
                                    color_list = NULL,
                                    peptide_annotation = NULL,
                                    feature_id_col = 'peptide_group_label',
                                    factors_of_feature_ann = c('KEGG_pathway',
                                                               'evolutionary_distance'),
                                    color_list_features = NULL,
                                    filename = NULL, width = 7, height = 7, 
                                    units = c('cm','in','mm'), 
                                    plot_title = NULL,
                                    ...){
  
  df_long = matrix_to_long(data_matrix, sample_id_col = sample_id_col)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     merge = FALSE)
  data_matrix = long_to_matrix(df_long, sample_id_col = sample_id_col)
  rm(df_long)
  
  #infer the color scheme for sample annotation (cols)
  if(is.null(color_list) && !is.null(sample_annotation)){
    warning('color_list for samples (cols) not defined, inferring automatically.
            Numeric/factor columns are guessed, for more controlled color 
            mapping use sample_annotation_to_colors()')
    color_list = sample_annotation_to_colors(sample_annotation = sample_annotation, 
                                             sample_id_col = sample_id_col, 
                                             factor_columns = factors_to_plot,
                                             numeric_columns = NULL,
                                             guess_factors = TRUE)
  }
  
  #infer the color scheme for feature annotation (rows)
  if(is.null(color_list_features) && !is.null(peptide_annotation)){
    warning('color_list_features for features (rows) not defined, inferring 
            automatically. Numeric/factor columns are guessed, for more 
            controlled color mapping use sample_annotation_to_colors()')
    color_list_features = sample_annotation_to_colors(peptide_annotation, 
                                                      sample_id_col = feature_id_col, 
                                                      factor_columns = factors_of_feature_ann,
                                                      numeric_columns = NULL,
                                                      guess_factors = TRUE)
  }
  
  p <- plot_heatmap_generic(data_matrix, 
                            column_annotation_df = sample_annotation,
                            row_annotation_df = peptide_annotation, 
                            fill_the_missing = fill_the_missing, 
                            col_ann_id_col = sample_id_col,
                            row_ann_id_col = feature_id_col,
                            columns_for_cols = factors_to_plot,
                            columns_for_rows = factors_of_feature_ann,
                            cluster_rows = cluster_cols, 
                            cluster_cols = cluster_cols,
                            annotation_color_cols = color_list,
                            annotation_color_rows = color_list_features,
                            heatmap_color = heatmap_color,
                            color_for_missing = color_for_missing,
                            filename = filename, width = width, height = width, 
                            units = units, 
                            plot_title = plot_title,
                            ...)
  return(p)
}

#' Plot the heatmap
#'
#' @inheritParams proBatch
#' @param data_matrix the matrix of data to be plotted
#' @param column_annotation_df data frame annotating columns of 
#' \code{data_matrix}
#' @param row_annotation_df data frame annotating rows of \code{data_matrix}
#' @param col_ann_id_col column of \code{column_annotation_df} whose values are
#' unique identifiers of columns in \code{data_matrix}
#' @param row_ann_id_col column of \code{row_annotation_df} whose values are 
#' unique identifiers of rows in \code{data_matrix}
#' @param columns_for_cols vector of factors (columns) of 
#' \code{column_annotation_df}
#' that will be mapped to color annotation of heatmap columns
#' @param columns_for_rows vector of factors (columns) of 
#' \code{row_annotation_df}
#' that will be mapped to color annotation of heatmap rows
#' @param cluster_rows boolean: whether the rows should be clustered
#' @param cluster_cols boolean: whether the rows should be clustered
#' @param annotation_color_cols list of color vectors for column annotation,
#' for each factor to be plotted; for factor-like variables a named vector 
#' (names should correspond to the levels of factors). Advisable to supply here
#' color list returned by \code{sample_annotation_to_colors}
#' @param annotation_color_rows list of color vectors for row annotation,
#' for each factor to be plotted; for factor-like variables a named vector 
#' (names should correspond to the levels of factors). Advisable to supply here
#' color list returned by \code{sample_annotation_to_colors}
#' @param fill_the_missing numeric value that the missing values are
#'   substituted with, or \code{NULL} if features with missing values are to be 
#'   excluded.
#' @param color_for_missing special color to make missing values. 
#' Usually black or white, depending on \code{heatmap_color}
#' @param heatmap_color vector of colors used in heatmap (typicall a gradient)
#' @param ... other parameters of \code{link[pheatmap]{pheatmap}}
#'
#' @return pheatmap-type object
#' @export
#' 
#' @examples 
#' 
#' p <- plot_heatmap_generic(log_transform_dm(example_proteome_matrix), 
#' column_annotation_df = example_sample_annotation,
#' columns_for_cols = c("MS_batch",  "digestion_batch", "Diet", 'DateTime'),
#' plot_title = 'test_heatmap',
#' show_rownames = FALSE, show_colnames = FALSE)
#' 
plot_heatmap_generic <- function(data_matrix, 
                                 column_annotation_df = NULL,
                                 row_annotation_df = NULL, 
                                 col_ann_id_col = 'FullRunName',
                                 row_ann_id_col = 'peptide_group_label',
                                 columns_for_cols = c('MS_batch','Diet', 
                                                      'DateTime','order'),
                                 columns_for_rows = c('KEGG_pathway',
                                                      'WGCNA_module',
                                                      'evolutionary_distance'),
                                 cluster_rows = FALSE, cluster_cols = TRUE,
                                 annotation_color_cols = NULL,
                                 annotation_color_rows = NULL,
                                 fill_the_missing = -1, 
                                 color_for_missing = 'black',
                                 heatmap_color = colorRampPalette(
                                   rev(brewer.pal(n = 7, 
                                                  name = "RdYlBu")))(100),
                                 filename = NULL, width = 7, height = 7, 
                                 units = c('cm','in','mm'), 
                                 plot_title = NULL,
                                 ...){
  
  #deal with the missing values
  warning_message <- 'Heatmap cannot operate with missing values in the matrix'
  data_matrix = handle_missing_values(data_matrix, warning_message, 
                                      fill_the_missing)
  if (is.null(fill_the_missing) & any(is.na(data_matrix)) &
      (!cluster_rows | !cluster_cols) ){
    message('With NAs removed, clustering of heatmap will work, 
              specify: cluster_rows = T, cluster_cols = T')
  }
  
  if(!is.null(fill_the_missing)){
    heatmap_color = c(color_for_missing, heatmap_color)
  }
  
  annotation_col = NA
  annotation_row = NA
  if(!is.null(column_annotation_df)){
    if(!is.null(columns_for_cols)){
      annotation_col = column_annotation_df %>% 
        select(one_of(col_ann_id_col, columns_for_cols))
    } else {
      annotation_col = column_annotation_df
    }
    annotation_col = annotation_col %>%
      mutate_if(is.POSIXct, as.numeric) %>%
      remove_rownames %>% 
      column_to_rownames(var=col_ann_id_col)  
  }
  
  if(!is.null(row_annotation_df)){
    if(!is.null(columns_for_rows)){
      annotation_row = row_annotation_df %>% 
        select(one_of(row_ann_id_col, columns_for_rows))
    } else {
      annotation_row = row_annotation_df
    }
    
    annotation_row = annotation_row %>%
      mutate_if(is.POSIXct, as.numeric) %>%
      remove_rownames %>% 
      column_to_rownames(var=row_ann_id_col)
  }
  
  if(is.null(column_annotation_df) && is.null(row_annotation_df)){
    warning("annotation_row and annotation_col are not specified for heatmap 
            (annotation of rows/cols such as sample annotation will not be plotted)")
  }
  
  if(!identical(annotation_col, NA) & (is.data.frame(annotation_col) | is.matrix(annotation_col))){
    if (!setequal(rownames(annotation_col), colnames(data_matrix))){
      warning('coloring by column annotation will not work: annotation rownames do not match data matrix column names')
    }
  }
  
  if(!identical(annotation_row, NA) & (is.data.frame(annotation_row) | is.matrix(annotation_row))){
    if (!setequal(rownames(annotation_row), rownames(data_matrix))){
      warning('coloring by row annotation will not work: annotation rownames do not match data matrix column names')
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
  
  if(is.null(annotation_color_cols) && is.null(annotation_color_rows)){
    annotation_color_list = NA
  } else {
    annotation_color_list = c(annotation_color_cols, annotation_color_rows)
  }
  
  p <- pheatmap(data_matrix, 
                cluster_rows = cluster_rows, cluster_cols = cluster_cols,
                color = heatmap_color,
                annotation_col = annotation_col, 
                annotation_row = annotation_row, 
                annotation_colors = annotation_color_list,
                filename = filename, width = width, height = height,
                main = plot_title, ...)
  return(p)
}

#' Calculate variance distribution by variable
#' 
#' @inheritParams proBatch
#' @param factors_for_PVCA vector of factors from \code{sample_annotation}, that
#'   are used in PVCA analysis
#' @param pca_threshold the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param variance_threshold the percentile value of weight each of the factors
#'   needs to explain (the rest will be lumped together)
#' @param fill_the_missing numeric value determining how  missing values 
#' should be substituted. If \code{NULL}, features with missing values are 
#' excluded.
#' @return data frame of weights of Principal Variance Components
#' @export
#' 
#' @examples
#' 
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_df <- calculate_PVCA(matrix_test, example_sample_annotation, 
#' factors_for_PVCA = c('MS_batch', 'digestion_batch',"Diet", "Sex", "Strain"),
#' pca_threshold = .6, variance_threshold = .01, fill_the_missing = -1)
calculate_PVCA <- function(data_matrix, sample_annotation, 
                           feature_id_col = 'peptide_group_label',
                           sample_id_col = 'FullRunName',
                           factors_for_PVCA = c('MS_batch', 'digestion_batch',
                                                "Diet", "Sex", "Strain"),
                           pca_threshold = .6, variance_threshold = .01,
                           fill_the_missing = -1) {
  
  df_long = matrix_to_long(data_matrix, sample_id_col = sample_id_col)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long,
                                     batch_col = NULL, order_col = NULL, 
                                     facet_col = NULL, merge = FALSE)
  data_matrix = long_to_matrix(df_long, sample_id_col = sample_id_col)
  
  sample_annotation = sample_annotation   %>% 
    select(one_of(c(sample_id_col, factors_for_PVCA))) %>%
    mutate_if(is.POSIXct, as.numeric) %>%
    as.data.frame()%>%
    column_to_rownames(var = sample_id_col)
  
  data_matrix = check_feature_id_col_in_dm(feature_id_col, data_matrix)
  
  warning_message <- 'PVCA cannot operate with missing values in the matrix'
  data_matrix = handle_missing_values(data_matrix, warning_message, 
                                      fill_the_missing)
  
  covrts.annodf = Biobase::AnnotatedDataFrame(data=sample_annotation)
  data_matrix = data_matrix[, rownames(sample_annotation)]
  expr_set = Biobase::ExpressionSet(assayData = data_matrix, 
                                    phenoData = covrts.annodf)
  pvcaAssess = pvcaBatchAssess(expr_set, factors_for_PVCA, 
                               threshold = pca_threshold)
  pvcaAssess_df = data.frame(weights = as.vector(pvcaAssess$dat),
                             label = pvcaAssess$label,
                             stringsAsFactors = FALSE)
  
  label_of_small = sprintf('Below %1.0f%%', 100*variance_threshold)
  if (sum(pvcaAssess_df$weights < variance_threshold) > 1){
    small_weights <- pvcaAssess_df$weights < variance_threshold
    pvca_res_small = sum(pvcaAssess_df$weights[small_weights])
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
#' @param technical_factors vector \code{sample_annotation} column names that 
#' are technical covariates
#' @param biological_factors vector \code{sample_annotation} column names, that
#'   are biologically meaningful covariates
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#' @param pca_threshold the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param variance_threshold the percentile value of weight each of the 
#'  covariates needs to explain (the rest will be lumped together)
#' @param fill_the_missing numeric value determining how  missing values 
#' should be substituted. If \code{NULL}, features with missing values are 
#' excluded.
#' If \code{NULL}, features with missing values are excluded.
#'
#' @return \code{ggplot} object with the plot
#' @export
#'
#' @examples 
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_plot <- plot_PVCA(matrix_test, example_sample_annotation, 
#' technical_factors = c('MS_batch', 'digestion_batch'),
#' biological_factors = c("Diet", "Sex", "Strain"))
#' 
#' \dontrun{
#' pvca_plot <- plot_PVCA(matrix_test, example_sample_annotation, 
#' technical_factors = c('MS_batch', 'digestion_batch'),
#' biological_factors = c("Diet", "Sex", "Strain"), 
#' filename = 'test_PVCA.png', width = 28, height = 22, units = 'cm')
#' }
#' 
#' @seealso \code{\link{sample_annotation_to_colors}}, 
#' \code{\link[ggplot2]{ggplot}}
plot_PVCA <- function(data_matrix, sample_annotation,
                      feature_id_col = 'peptide_group_label',
                      sample_id_col = 'FullRunName',
                      technical_factors = c('MS_batch', 'instrument'),
                      biological_factors = c('cell_line','drug_dose'),
                      fill_the_missing = -1,
                      pca_threshold = .6, variance_threshold = .01,
                      colors_for_bars = NULL,
                      filename = NULL, width = NA, height = NA, 
                      units = c('cm','in','mm'),
                      plot_title = NULL,
                      theme = 'classic',
                      base_size = 20){
  pvca_res = prepare_PVCA_df(data_matrix = data_matrix, 
                             sample_annotation = sample_annotation,
                             feature_id_col = feature_id_col,
                             sample_id_col = sample_id_col,
                             technical_factors = technical_factors,
                             biological_factors = biological_factors,
                             fill_the_missing = fill_the_missing,
                             pca_threshold = pca_threshold, 
                             variance_threshold = variance_threshold)
  
  gg = plot_PVCA.df(pvca_res = pvca_res, colors_for_bars = colors_for_bars, 
                    filename = filename, width = width, height = height, units = units, 
                    plot_title = plot_title,
                    theme = theme, base_size = base_size)
  return(gg)
}

#' prepare the weights of Principal Variance Components 
#'
#' @inheritParams proBatch
#' @param technical_factors vector \code{sample_annotation} column names that 
#' are technical covariates
#' @param biological_factors vector \code{sample_annotation} column names, that
#'   are biologically meaningful covariates
#' @param pca_threshold the percentile value of the minimum amount of the
#'   variabilities that the selected principal components need to explain
#' @param variance_threshold the percentile value of weight each of the 
#'  covariates needs to explain (the rest will be lumped together)
#' @param fill_the_missing numeric value determining how  missing values 
#' should be substituted. If \code{NULL}, features with missing values are 
#' excluded.
#' If \code{NULL}, features with missing values are excluded.
#'
#' @return data frame with weights and factors, combined in a way ready for plotting
#' @export
#'
#' @examples
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_df_res <- prepare_PVCA_df(matrix_test, example_sample_annotation, 
#' technical_factors = c('MS_batch', 'digestion_batch'),
#' biological_factors = c("Diet", "Sex", "Strain"), 
#' pca_threshold = .6, variance_threshold = .01, fill_the_missing = -1)
prepare_PVCA_df <- function(data_matrix, sample_annotation,
                            feature_id_col = 'peptide_group_label',
                            sample_id_col = 'FullRunName',
                            technical_factors = c('MS_batch', 'instrument'),
                            biological_factors = c('cell_line','drug_dose'),
                            fill_the_missing = -1,
                            pca_threshold = .6, variance_threshold = .01){
  
  factors_for_PVCA = c(technical_factors, biological_factors)
  
  pvca_res = calculate_PVCA(data_matrix, sample_annotation, 
                            feature_id_col = feature_id_col,
                            sample_id_col = sample_id_col,
                            factors_for_PVCA = factors_for_PVCA,
                            pca_threshold = pca_threshold, 
                            variance_threshold = variance_threshold,
                            fill_the_missing = fill_the_missing)
  
  tech_interactions = expand.grid(technical_factors, 
                                  technical_factors) %>%
    mutate(tech_interactions = paste(Var1, Var2, sep = ':')) %>%
    pull(tech_interactions)
  biol_interactions = expand.grid(biological_factors, 
                                  biological_factors) %>%
    mutate(biol_interactions = paste(Var1, Var2, sep = ':')) %>%
    pull(biol_interactions)
  
  label_of_small = sprintf('Below %1.0f%%', 100*variance_threshold)
  technical_factors = c(technical_factors, tech_interactions)
  biological_factors = c(biological_factors, biol_interactions)
  pvca_res = pvca_res %>% mutate(category = ifelse(label %in% technical_factors, 
                                                   'technical',
                                                   ifelse(label %in% biological_factors, 
                                                          'biological',
                                                          ifelse(label %in% c(label_of_small, 'resid'), 
                                                                 'residual', 'biol:techn'))))
  
  pvca_res = pvca_res %>%
    arrange(desc(weights)) %>%
    arrange(label == label_of_small) %>%
    arrange(label == 'resid')
  return(pvca_res)
}

#' plot PVCA, when the analysis is completed
#'
#' @inheritParams proBatch
#' @param pvca_res data frame of weights of Principal Variance Components, result
#' of \code{calculate_PVCA}
#' @param colors_for_bars four-item color vector, specifying colors for the
#'   following categories: c('residual', 'biological', 'biol:techn',
#'   'technical')
#'   
#' @return \code{ggplot} object with bars as weights, colored by bio/tech factors
#' @export
#'
#' @examples
#' matrix_test <- example_proteome_matrix[1:150, ]
#' pvca_df_res <- prepare_PVCA_df(matrix_test, example_sample_annotation, 
#' technical_factors = c('MS_batch', 'digestion_batch'),
#' biological_factors = c("Diet", "Sex", "Strain"), 
#' pca_threshold = .6, variance_threshold = .01, fill_the_missing = -1)
#' colors_for_bars = c('grey', 'green','blue','red')
#' names(colors_for_bars) = c('residual', 'biological','biol:techn','technical')
#' 
#' pvca_plot <- plot_PVCA.df(pvca_df_res, colors_for_bars)
plot_PVCA.df <- function(pvca_res,
                      colors_for_bars = NULL,
                      filename = NULL, width = NA, height = NA, 
                      units = c('cm','in','mm'),
                      plot_title = NULL,
                      theme = 'classic',
                      base_size = 20){
  pvca_res = pvca_res %>%
    mutate(label = factor(label, levels=label))
  
  y_title = 'Weighted average proportion variance'
  gg  = ggplot(pvca_res, aes(x = label, y = weights, fill = category))+
    geom_bar(stat = 'identity', color = 'black')+
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
    gg = gg + theme_classic(base_size = base_size)
  }else{
    message("plotting with default ggplot theme, only theme = 'classic' 
            implemented")
  }
  
  save_ggplot(filename, units, width, height, gg)
  
  gg = gg +
    theme(axis.title.x = NULL, 
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))+
    xlab(NULL)+
    theme(text = element_text(size=15))+
    guides(fill=guide_legend(override.aes=list(color=NA), title=NULL))
  
  return(gg)
}

#' plot PCA plot
#'
#' @inheritParams proBatch
#' @param color_by column name (as in \code{sample_annotation}) to color by
#' @param PC_to_plot principal component numbers for x and y axis
#' @param fill_the_missing numeric value determining how  missing values 
#' should be substituted. If \code{NULL}, features with missing values are 
#' excluded.
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
#' color_by = 'DateTime', plot_title = "PCA colored by DateTime")
#' 
#' color_list <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch', 'digestion_batch'),
#' numeric_columns = c('DateTime','order'))
#' pca_plot <- plot_PCA(example_proteome_matrix, example_sample_annotation, 
#' color_by = 'DateTime', color_scheme = color_list[['DateTime']])
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
                     sample_id_col = 'FullRunName',
                     color_by = 'MS_batch',
                     PC_to_plot = c(1,2), fill_the_missing = -1,
                     color_scheme = 'brewer',
                     filename = NULL, width = NA, height = NA, 
                     units = c('cm','in','mm'),
                     plot_title = NULL,
                     theme = 'classic',
                     base_size = 20){
  
  df_long = matrix_to_long(data_matrix, sample_id_col = sample_id_col)
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col = color_by, order_col = NULL, 
                                     facet_col = NULL, merge = FALSE)
  data_matrix = long_to_matrix(df_long, sample_id_col = sample_id_col)
  
  data_matrix = check_feature_id_col_in_dm(feature_id_col, data_matrix)
  
  warning_message <- 'PCA cannot operate with missing values in the matrix'
  data_matrix = handle_missing_values(data_matrix, warning_message, 
                                      fill_the_missing)
  
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
    gg = gg + theme_classic(base_size = base_size)
  }else{
    message("plotting with default ggplot theme, only theme = 'classic' 
            implemented")
  }
  
  save_ggplot(filename, units, width, height, gg)
  
  return(gg)
}


