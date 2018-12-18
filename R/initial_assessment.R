#' Plot per-sample mean or boxplot (showing median and quantiles) vs order (if the real
#' running order available)
#' @details functions for quick visual assessment of trends associated, overall
#'   or specific covariate-associated (see `batch_col` and `facet_col`)
#' @param data_matrix features (in rows) vs samples (in columns) matrix, 
#' with feature IDs in rownames and file/sample names as colnames. in most function,
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
#' @param vline_color color of vertical line, typically to denote batches
#' @param ylimits range of y-axis to plot feature-level trends 
#' @return ggplot2 class object. Thus, all aesthetics can be overriden
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @name plot_sample_mean_or_boxplot
#'
#' @export
#'
#' @examples \dontrun{plot_sample_mean(log_transformed_matrix, example_sample_annotation, 
#' order_col = 'order', batch_col = "MS_batch", ylimits = c(12, 16))}
#' 
plot_sample_mean <- function(data_matrix, sample_annotation = NULL,
                             sample_id_col = 'FullRunName',
                             order_col = 'order',
                             batch_col = "MS_batch",
                             facet_col = NULL,
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
#' @examples \dontrun{plot_boxplot(log_transformed_long, example_sample_annotation, 
#' batch_col = "MS_batch"}
#' 
plot_boxplot <- function(df_long, sample_annotation = NULL,
                         sample_id_col = 'FullRunName',
                         measure_col = 'Intensity',
                         order_col = 'order',
                         batch_col = 'MS_batch',
                         facet_col = NULL,
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
