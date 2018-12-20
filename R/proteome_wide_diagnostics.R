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
#' @param label_samples if \code{TRUE} sample IDs (column names of 
#' \code{data_matrix}) will be printed
#' @param label_font size of the font. Is active if \code{label_samples} is 
#' \code{TRUE}, ignored otherwise
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param ... other parameters of `plotDendroAndColors` from `WGCNA` package
#'
#' @return No return
#' @examples
#' color_scheme <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch','EarTag', "Strain", "Diet", "digestion_batch", "Sex"),
#' not_factor_columns = 'DateTime',
#' numeric_columns = c('order'))
#' 
#' color_annotation <- color_scheme$color_df
#' 
#' plot_hierarchical_clustering(example_proteome_matrix, color_annotation,  
#' distance = "euclidean", agglomeration = 'complete',
#' label_samples = FALSE)
#' 
#' @export
#'
#' @seealso \code{\link[stats]{hclust}}, \code{\link{sample_annotation_to_colors}},
#'   \code{\link[WGCNA]{plotDendroAndColors}}
plot_hierarchical_clustering  <- function(data_matrix, color_df,
                                          distance = "euclidean",
                                          agglomeration = 'complete',
                                          label_samples = TRUE, label_font = .2,
                                          plot_title = NULL,
                                          ...){
    dist_matrix = dist(t(as.matrix(data_matrix)), method = distance)
    hierarchical_clust = hclust(dist_matrix, method = agglomeration)
    if (label_samples){
        if (ncol(data_matrix) > 80){
            warning('Too many samples, adjust the font with `label_font` argument or
              remove labels by setting `label_samples = FALSE` in function call')
        }
        WGCNA::plotDendroAndColors(hierarchical_clust, color_df, rowTextAlignment = 'left',
                                   main = plot_title,
                                   hang = -0.1, addGuide = TRUE, dendroLabels = NULL,
                                   cex.dendroLabels = label_font, ...)

    } else{
        WGCNA::plotDendroAndColors(hierarchical_clust, color_df, rowTextAlignment = 'left',
                                   main = plot_title,
                                   hang = -0.1, addGuide = TRUE, dendroLabels = FALSE, ...)
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
#' }; each column of sample annotation will get it's own row. 
#' If `cluster_cols = T` this will indicate,
#' whether sample proximity is driven by one of 
#' biolical or technical factors
#' @param sample_id_col name of the column in 
#' sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param fill_the_missing boolean value determining if 
#' missing values should be
#'   substituted with -1 (and colored with black)
#' @param cluster_rows boolean value determining if rows 
#' should be clustered
#' @param cluster_cols boolean value determining if columns 
#' should be clustered
#' @param sample_annotation_col biological or technical 
#' factors to be annotated in heatmap columns
#' @param sample_annotation_row biological or technical 
#' factors to be annotated in heatmap rows 
#' @param annotation_color_list list specifying colors 
#' for columns (samples). Best created by `sample_annotation_to_colors`
#' @param heatmap_color vector of colors used in heatmap (typicall a gradient)
#' @param color_for_missing special color to make missing values. 
#' Usually black or white, depending on `heatmap_color`
#' @param filename filepath where to save the image
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param ... other parameters of \code{link[pheatmap]{pheatmap}}
#' 
#' @return object returned by \code{link[pheatmap]{pheatmap}}
#' @export
#' 
#' @examples 
#' color_scheme <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch','EarTag', "Strain", 
#' "Diet", "digestion_batch", "Sex"),
#' not_factor_columns = 'DateTime',
#' numeric_columns = c('order'))
#' 
#' plot_heatmap(example_proteome_matrix, 
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
                         fill_the_missing = TRUE, cluster_rows = TRUE, cluster_cols = FALSE,
                         annotation_color_list = NA,
                         heatmap_color = colorRampPalette(
                           rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(100),
                         color_for_missing = 'black',
                         filename = NA, plot_title = NA,
                         ...){
    if(fill_the_missing) {
        data_matrix[is.na(data_matrix)] = 0
        warning('substituting missing values with 0, 
                this might affect the clustering!')
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

    covrts.annodf = Biobase::AnnotatedDataFrame(data=sample_annotation)
    expr_set = Biobase::ExpressionSet(data_matrix[,rownames(sample_annotation)], covrts.annodf)
    pvcaAssess = pvca::pvcaBatchAssess (expr_set, factors_for_PVCA, threshold = threshold_pca)
    pvcaAssess_df = data.frame(weights = as.vector(pvcaAssess$dat),
                               label = pvcaAssess$label,
                               stringsAsFactors = FALSE)

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
#' @param fill_the_missing numeric value that the missing values are
#'   substituted with
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see
#'   examples)
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#'
#' @return list of two items: plot =gg, df = pvca_res
#' @export
#'
#' @examples \dontrun{plot_PVCA(example_proteome_marix, example_sample_annotation, 
#' technical_covariates = c('MS_batch', 'digestion_batch'),
#' biological_covariates = c("Diet", "Sex", "Strain"))}
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

    if(setequal(unique(sample_annotation[[sample_id_col]]), 
                unique(colnames(data_matrix))) == FALSE){
        warning('Sample IDs in sample annotation not 
                consistent with samples in input data.')}
    
    sample_names = sample_annotation[[sample_id_col]]
    sample_annotation = sample_annotation %>% 
      select(one_of(factors_for_PVCA)) %>%
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
                        assign values of this column to rowname and remove from the data frame!', 
                                feature_id_col))
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
    if(theme == 'classic'){
        gg = gg +theme_classic()
    }
    if(is.null(colors_for_bars)){
        colors_for_bars = c('grey', wesanderson::wes_palettes$Rushmore[3:5])
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
#' plot_PCA(example_proteome_matrix, example_sample_annotation, 
#' color_by = 'MS_batch', plot_title = "PCA colored by MS batch")
#' 
#' @seealso \code{\link[ggfortify]{autoplot.pca_common}}, 
#' \code{\link[ggplot2]{ggplot}}
plot_PCA <- function(data_matrix, sample_annotation,
                     feature_id_col = 'peptide_group_label',
                     color_by = 'MS_batch',
                     PC_to_plot = c(1,2), fill_the_missing = 0,
                     colors_for_factor = NULL,
                     theme = 'classic',
                     plot_title = NULL){

    if(!is.null(feature_id_col)){
        if(feature_id_col %in% colnames(data_matrix)){
            if(is.data.frame(data_matrix)){
                warning(sprintf('feature_id_col with name %s in data matrix instead of rownames,
                        this might cause errors in other diagnostic functions,
                        assign values of this column to rowname and remove from the data frame!', 
                                feature_id_col))
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


