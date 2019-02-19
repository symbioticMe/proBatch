#' Plot per-sample mean or boxplots (showing median and quantiles). In ordered samples,
#' e.g. consecutive MS runs, order-associated effects are visualised.
#' @details functions for quick visual assessment of trends associated, overall
#'   or specific covariate-associated (see \code{batch_col} and \code{facet_col})
#' @param data_matrix features (in rows) vs samples (in columns) matrix, 
#' with feature IDs in rownames and file/sample names as colnames. in most function,
#' @param df_long data frame where each row is a single feature in a single
#'   sample, thus it has minimally, \code{sample_id_col}, \code{feature_id_col} and
#'   \code{measure_col}, but usually also \code{m_score} (in OpenSWATH output result
#'   file)
#' @param sample_annotation data matrix with 1) \code{sample_id_col} (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param measure_col if \code{df_long} is among the parameters, it is the column
#'   with expression/abundance/intensity, otherwise, it is used internally for
#'   consistency
#' @param batch_col column in \code{sample_annotation} that should be used for
#'   batch comparison. Can be `NULL` if only boxplot/mean comparison, without coloring by batches, is required.
#' @param order_col column where running order is specified.
#' @param color_by_batch should the each batch be represented with its own
#'   color?
#' @param color_scheme named vector, names corresponding to unique batch values
#'   as specified in \code{sample_annotation}
#' @param facet_col column  in \code{sample_annotation} with a batch factor to separate 
#' plots into facets; usually 2nd to \code{batch_col}. Most meaningful for multi-instrument 
#' MS experiments (where each instrument has its own order-associated effects) 
#' or simultaneous examination of two batch factors (e.g. preparation day and measurement day)
#' @param theme ggplot theme, by default  \code{classic}. Can be easily overriden (see
#'   examples)
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param order_per_facet if order is defined ignoring facets (usually
#'   instrument), re-define order per-batch
#' @param vline_color color of vertical lines, typically denoting 
#'  different MS batches in ordered runs; should be \code{NULL} for experiments without intrinsic order
#' @param ylimits range of y-axis to plot feature-level trends 
#' @return ggplot2 class object. Thus, all aesthetics can be overriden
#'
#' @seealso \code{\link[ggplot2]{ggplot}}
#' @name plot_sample_mean_or_boxplot
#'
#' @export
#'
#' @examples 
#' plot_sample_mean(example_proteome_matrix, example_sample_annotation, 
#' order_col = 'order', batch_col = "MS_batch")
#' 
plot_sample_mean <- function(data_matrix, sample_annotation = NULL,
                             sample_id_col = 'FullRunName',
                             batch_col = "MS_batch",
                             color_by_batch = FALSE, color_scheme = 'brewer',
                             order_col = 'order',
                             vline_color = 'grey',
                             facet_col = NULL,
                             plot_title = NULL,
                             theme = 'classic',
                             ylimits = NULL){
  
  #Create a data frame with sample averages
  sample_average = colMeans(data_matrix, na.rm = TRUE)
  #names(sample_average) = colnames(data_matrix)
  
  df_ave = data.frame(Mean_Intensity = sample_average,
                      sample_id_col = colnames(data_matrix), 
                      stringsAsFactors = FALSE)
  
  #Assign a column for sample_id
  if(is.null(sample_annotation) | !(sample_id_col %in% names(sample_annotation))
     | is.null(sample_id_col)){
    sample_id_col = 'data_matrix_colnames'
  }
  names(df_ave)[names(df_ave) == "sample_id_col"] <- sample_id_col
  
  #Check the consistency of sample annotation sample IDs and measurement table sample IDs
  df_ave = check_sample_consistency(sample_annotation, sample_id_col, df_ave)
  
  #Ensure that batch-coloring-related arguments are defined properly
  if(!is.null(batch_col)){
    if(!(batch_col %in% names(df_ave))){
      stop('batches cannot be colored as the batch column or sample ID column
             is not defined, check sample_annotation and data matrix')
    }
    #For proper plotting, batch column has to be a factor
    df_ave[, batch_col] <- as.factor(df_ave[, batch_col])
  } else {
    if (color_by_batch){
      warning('batches cannot be colored as the batch column is defined as NULL, continuing without colors')
      color_by_batch = FALSE
    }
  }
  
  #For order definition and subsequent faceting, facet column has to be in the data frame
  if(!is.null(facet_col)){
    if ( !(facet_col %in% names(df_ave))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present in the data,
                 check sample annotation data frame', facet_col))
    }
  }
  
  #Defining sample order for plotting
  sample_order = define_sample_order(order_col, sample_annotation, facet_col, batch_col, df_ave, 
                      sample_id_col, color_by_batch)
  order_col = sample_order$order_col
  df_ave = sample_order$df_long
  
  #Main plotting of intensity means:
  gg = ggplot(df_ave, aes_string(x = order_col, y = 'Mean_Intensity'))+
    geom_point()
  
  #add colors
  gg = color_points_by_batch(color_by_batch, batch_col, gg, color_scheme, sample_annotation)
  
  #add vertical lines, if required (for order-related effects)
  gg = add_vertical_batch_borders(order_col, sample_id_col, batch_col, vline_color, facet_col, df_ave, gg)
  
  #Plot each "facet factor" in it's own subplot
  if(!is.null(facet_col)){
    gg = gg + facet_wrap(as.formula(paste("~", facet_col)),
                         dir = 'v', scales = "free_x")
  }
  
  #Change the theme
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  
  #Add the title
  if(!is.null(plot_title)) {
    gg = gg + ggtitle(plot_title)+
    theme(plot.title = element_text(face = 'bold',hjust = .5))
  }
  
  #Rotate x axis tick labels if the filenames, not numeric order, is displayed
  if (!is.numeric(df_ave[[order_col]])){
    if(is.character(df_ave[[order_col]])){
      df_ave[[order_col]] = factor(df_ave[[order_col]],
                                   levels = unique(df_ave[[order_col]]))
    }
    gg = gg +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  }
  
  #Move the legend to the upper part of the plot to save the horizontal space
  if (length(unique(df_ave[[order_col]])) > 30){
    gg = gg + theme(legend.position="top")
  }
  
  #Change the limits of vertical axes
  if(!is.null(ylimits)){
    gg = gg +
      ylim(ylimits)
  }
  
  return(gg)
}


#' @name plot_sample_mean_or_boxplot
#' 
#' @export
#'
#' @examples
#' plot_boxplot(example_proteome, example_sample_annotation, 
#' batch_col = "MS_batch")
#' 
plot_boxplot <- function(df_long, sample_annotation = NULL,
                         sample_id_col = 'FullRunName',
                         measure_col = 'Intensity',
                         batch_col = 'MS_batch',
                         color_by_batch = TRUE, color_scheme = 'brewer',
                         order_col = 'order',
                         facet_col = NULL,
                         plot_title = NULL, theme = 'classic',
                         ylimits = NULL){
  
  #Check the consistency of sample annotation sample IDs and measurement table sample IDs
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long)
  
  #Ensure that batch-coloring-related arguments are defined properly
  if(!is.null(batch_col)){
    if(!(batch_col %in% names(df_long))){
      stop('batches cannot be colored as the batch column or sample ID column
             is not defined, check sample_annotation and data matrix')
    }
    #For proper plotting, batch column has to be a factor
    df_long[, batch_col] <- as.factor(df_long[, batch_col])
  } else {
    if (color_by_batch){
      warning('batches cannot be colored as the batch column is defined as NULL, continuing without colors')
      color_by_batch = FALSE
    }
  }
  
  #For order definition and subsequent faceting, facet column has to be in the data frame
  if(!is.null(facet_col)){
    if ( !(facet_col %in% names(df_long))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present in the data,
                 check sample annotation data frame', facet_col))
    }
  }
  
  #Defining sample order for plotting
  sample_order = define_sample_order(order_col, sample_annotation, facet_col, batch_col, df_long, 
                                     sample_id_col, color_by_batch)
  order_col = sample_order$order_col
  df_long = sample_order$df_long
  
  #Main plotting of intensity distribution boxplots
  gg = ggplot(df_long, aes_string(x = order_col, y = measure_col,
                                  group = order_col))+
    geom_boxplot(outlier.size = .15)
  
  #Define the color scheme, add colors
  gg = color_fill_boxes_by_batch(color_by_batch, batch_col, gg, color_scheme, df_long)
  
  #Plot each "facet factor" in it's own subplot
  if(!is.null(facet_col)){
    gg = gg + facet_wrap(as.formula(paste("~", facet_col)),
                         dir = 'v', scales = "free_x")
  }
  
  #Change the plot theme
  if(theme == 'classic'){
    gg = gg + theme_classic()
  }
  
  #Add the title
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)+
      theme(plot.title = element_text(hjust = .5, face = 'bold', size = 16))
  }
  
  #Rotate x axis tick labels if the filenames, not numeric order, is displayed
  if (!is.numeric(df_long[[order_col]])){
    if(is.character(df_long[[order_col]])){
      df_long[[order_col]] = factor(df_long[[order_col]],
                                    levels = unique(df_long[[order_col]]))
    }
    gg = gg +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  }
  
  #Move the legend to the upper part of the plot to save the horizontal space
  if (length(unique(df_long[[order_col]])) > 30){
    gg = gg + theme(legend.position="top")
  }
  
  #Change the limits of vertical axes
  if(!is.null(ylimits)){
    gg = gg +
      ylim(ylimits)
  }
  
  return(gg)
}
