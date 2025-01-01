#' @title Plot per-sample mean or boxplots for initial assessment
#' @description Plot per-sample mean or boxplots (showing median and quantiles). In ordered samples,
#' e.g. consecutive MS runs, order-associated effects are visualised.
#' @details functions for quick visual assessment of trends associated, overall
#'   or specific covariate-associated (see \code{batch_col} and \code{facet_col})
#' @inheritParams proBatch
#' @param color_scheme named vector, names corresponding to unique batch values of 
#'  \code{batch_col} in \code{sample_annotation}. Best created with \link{sample_annotation_to_colors}
#' @param vline_color color of vertical lines, typically denoting 
#'  different MS batches in ordered runs; should be \code{NULL} for experiments without intrinsic order
#' @param ylimits range of y-axis to compare two plots side by side, if required.
#' @param outliers keep (default) or remove the boxplot outliers
#' 
#' @return ggplot2 class object. Thus, all aesthetics can be overridden
#'
#' @seealso \code{\link[ggplot2]{ggplot}}, \link{date_to_sample_order}
#' @name plot_sample_mean_or_boxplot
#'
#' @export
#'
#' @examples 
#' mean_plot <- plot_sample_mean(example_proteome_matrix, example_sample_annotation, 
#' order_col = 'order', batch_col = "MS_batch")
#' 
#' color_list <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch'),
#' numeric_columns = c('DateTime', 'order'))
#' plot_sample_mean(example_proteome_matrix, example_sample_annotation, 
#' order_col = 'order', batch_col = "MS_batch", color_by_batch = TRUE, 
#' color_scheme = color_list[["MS_batch"]])
#' 
#' \dontrun{
#' mean_plot <- plot_sample_mean(example_proteome_matrix, 
#'                               example_sample_annotation, 
#'                               order_col = 'order', batch_col = "MS_batch", 
#'                               filename = 'test_meanplot.png', 
#'                               width = 28, height = 18, units = 'cm')
#' }
#' 
plot_sample_mean <- function(data_matrix, sample_annotation = NULL,
                             sample_id_col = 'FullRunName',
                             batch_col = "MS_batch",
                             color_by_batch = FALSE, color_scheme = 'brewer',
                             order_col = 'order',
                             vline_color = 'grey',
                             facet_col = NULL,
                             filename = NULL, width = NA, height = NA, 
                             units = c('cm','in','mm'),
                             plot_title = NULL,
                             theme = 'classic',
                             base_size = 20,
                             ylimits = NULL){
  
  #Create a data frame with sample averages
  sample_average = colMeans(data_matrix, na.rm = TRUE)
  #names(sample_average) = colnames(data_matrix)
  
  df_ave = data.frame(Mean_Intensity = sample_average,
                      sample_id_col = colnames(data_matrix))
  names(df_ave)[names(df_ave) == "sample_id_col"] <- sample_id_col
  
  #Check the consistency of sample ann. sample IDs and measur. table sample IDs
  df_ave = check_sample_consistency(sample_annotation, sample_id_col, df_ave, 
                                    batch_col, order_col, facet_col)
  
  #Ensure that batch-coloring-related arguments are defined properly
  if(!is.null(batch_col)){
    if(!(batch_col %in% names(df_ave))){
      stop('batches cannot be colored as the batch column or sample ID column
             is not defined, check sample_annotation and data matrix')
    }
  } else {
    if (color_by_batch){
      warning('batches cannot be colored as the batch column is defined as NULL,
              continuing without colors')
      color_by_batch = FALSE
    }
  }
  
  #For order definition and subsequent faceting, facet column has to be in the 
  #data frame
  if(!is.null(facet_col)){
    if ( !(facet_col %in% names(df_ave))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present
                    in the data, check sample annotation data frame', 
                   facet_col))
    }
  }
  
  #Defining sample order for plotting
  sample_order = define_sample_order(order_col, sample_annotation, facet_col, 
                                     batch_col, df_ave, 
                      sample_id_col, color_by_batch)
  order_col = sample_order$order_col
  df_ave = sample_order$df_long
  
  #Main plotting of intensity means:
  gg = ggplot(df_ave, aes(x = !!sym(order_col), y = .data$Mean_Intensity))+
    geom_point()
  
  #add colors
  gg = color_by_factor(color_by_batch = color_by_batch, 
                       batch_col = batch_col, gg = gg, 
                       color_scheme = color_scheme, 
                       sample_annotation = df_ave,
                       fill_or_color = 'color')
  
  #add vertical lines, if required (for order-related effects)
  if (!is.null(batch_col)){
    batch_vector <- sample_annotation[[batch_col]]
    is_factor = is_batch_factor(batch_vector, color_scheme)
  }
  
  if (!is.null(batch_col) && is_factor){
    if (!is.null(sample_annotation)){
      gg = add_vertical_batch_borders(order_col, sample_id_col, batch_col, 
                                      vline_color, 
                                      facet_col, sample_annotation, gg)
    } else {
      gg = add_vertical_batch_borders(order_col, sample_id_col, batch_col, 
                                      vline_color, 
                                      facet_col, df_ave, gg)
    }
  }
  
  #Plot each "facet factor" in it's own subplot
  if(!is.null(facet_col)){
    gg = gg + facet_wrap(as.formula(paste("~", facet_col)),
                         dir = 'v', scales = "free_x")
  }
  
  #Add the title
  if(!is.null(plot_title)) {
    gg = gg + ggtitle(plot_title)+
      theme(plot.title = element_text(face = 'bold',hjust = .5))
  }
  
  #Change the theme
  if(!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic(base_size = base_size)
  } else{
    message("plotting with default ggplot theme, only theme = 'classic' 
            implemented")
  }
  
  #Change the limits of vertical axes
  if(!is.null(ylimits)){
    gg = gg +
      ylim(ylimits)
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
  
  if (length(unique(df_ave[[order_col]])) > 30 && color_by_batch && is_factor){
    gg = gg + theme(legend.position="top")
  }
  
  #save the plot
  save_ggplot(filename, units, width, height, gg)
  
  return(gg)
}


#' @name plot_sample_mean_or_boxplot
#' 
#' @export
#'
#' @examples
#' boxplot <- plot_boxplot(log_transform_df(example_proteome), 
#' sample_annotation = example_sample_annotation, 
#' batch_col = "MS_batch")
#' 
#' color_list <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch'),
#' numeric_columns = c('DateTime', 'order'))
#' plot_boxplot(log_transform_df(example_proteome), 
#' sample_annotation = example_sample_annotation, 
#' batch_col = "MS_batch", color_scheme = color_list[["MS_batch"]])

#' 
#' \dontrun{
#' boxplot <- plot_boxplot(log_transform_df(example_proteome), 
#' sample_annotation = example_sample_annotation, 
#' batch_col = "MS_batch", filename = 'test_boxplot.png', 
#' width = 14, height = 9, units = 'in')
#' }
#' 
plot_boxplot <- function(df_long, sample_annotation = NULL,
                         sample_id_col = 'FullRunName',
                         measure_col = 'Intensity',
                         batch_col = 'MS_batch',
                         color_by_batch = TRUE, color_scheme = 'brewer',
                         order_col = 'order',
                         facet_col = NULL,
                         filename = NULL, width = NA, height = NA, 
                         units = c('cm','in','mm'),
                         plot_title = NULL, theme = 'classic',
                         base_size = 20,
                         ylimits = NULL, outliers = TRUE){
  
  #Check the consistency of sample ann. sample IDs and measur. table sample IDs
  df_long = check_sample_consistency(sample_annotation, sample_id_col, df_long, 
                                     batch_col, order_col, facet_col)
  
  #Ensure that batch-coloring-related arguments are defined properly
  if(!is.null(batch_col)){
    if(!(batch_col %in% names(df_long))){
      stop('batches cannot be colored as the batch column or sample ID column
             is not defined, check sample_annotation and data matrix')
    }
  } else {
    if (color_by_batch){
      warning('batches cannot be colored as the batch column is defined as NULL, 
              continuing without colors')
      color_by_batch = FALSE
    }
  }
  
  #For order definition and subsequent faceting, facet column has to be in the 
  #data frame
  if(!is.null(facet_col)){
    if ( !(facet_col %in% names(df_long))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present
                    in the data, check sample annotation data frame', 
                   facet_col))
    }
  }
  
  #Defining sample order for plotting (even if order_col NULL, 
  #it will re-arrange df_long levels as required for plotting)
  sample_order = define_sample_order(order_col, sample_annotation, facet_col, 
                                     batch_col, df_long, 
                                     sample_id_col, color_by_batch)
  order_col = sample_order$order_col
  df_long = sample_order$df_long
  
  #Main plotting of intensity distribution boxplots
  gg = ggplot(df_long, aes(x = !!sym(order_col), y = !!sym(measure_col),
                           group = !!sym(order_col)))
  if (outliers){
    gg = gg + geom_boxplot(outlier.size = .15)
  } else {
    gg = gg + geom_boxplot(outlier.shape = NA)
  }
    
  
  #Define the color scheme, add colors
  gg = color_by_factor(color_by_batch = color_by_batch, 
                       batch_col = batch_col, gg = gg, 
                       color_scheme = color_scheme, 
                       sample_annotation = df_long,
                       fill_or_color = 'fill')

  #Plot each "facet factor" in it's own subplot
  if(!is.null(facet_col)){
    gg = gg + facet_wrap(as.formula(paste("~", facet_col)),
                         dir = 'v', scales = "free_x")
  }
  

  #Add the title
  if (!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)+
      theme(plot.title = element_text(hjust = .5, face = 'bold', size = 16))
  }
  
  #Change the plot theme
  if(!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic(base_size = base_size)
  } else{
    message("plotting with default ggplot theme, only theme = 'classic' 
            implemented")
  }
  
  #Change the limits of vertical axes
  if(!is.null(ylimits)){
    gg = gg +
      ylim(ylimits)
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
  
  if (!is.null(batch_col)){
    batch_vector <- sample_annotation[[batch_col]]
    is_factor = is_batch_factor(batch_vector, color_scheme)
  }
  
  #Move the legend to the upper part of the plot to save the horizontal space
  if (length(unique(df_long[[order_col]])) > 30  && 
      color_by_batch && is_factor){
    gg = gg + theme(legend.position="top")
  }
  
  #save the plot
  save_ggplot(filename, units, width, height, gg)
  
  return(gg)
}
