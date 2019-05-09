#' Plot peptide measurements
#'
#' Creates a peptide facetted ggplot2 plot of the value in \code{measure_col}
#' vs \code{order_col} (if `NULL`, x-axis is simply a sample name order). 
#' Additionally, the resulting plot can also be facetted by a batch factor, 
#' e.g. an MS instrument.
#'
#' @inheritParams proBatch
#' @param feature_name name of the selected feature (e.g. peptide) for diagnostic profiling
#' @param geom whether to show the feature as points and/or connect by lines (accepted 
#' values are: 1. \code{point}, \code{line} and \code{c('point', 'line')})
#' @param vline_color color of vertical lines, typically separating 
#'  different MS batches in ordered runs; 
#'  should be `NULL` for experiments without intrinsic order.
#' @param ylimits range of y-axis to plot feature-level trends 
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{feature_name} and (optionally) by \code{batch_col}
#' @examples 
#' single_feature_plot <- plot_single_feature(feature_name = "46213_NVGVSFYADKPEVTQEQK_2", 
#' df_long = example_proteome, example_sample_annotation, 
#' color_by_col = NULL)
#' 
#' #to examine features that have missing values specific per batch, this can be used:
#' plot_single_feature(feature_name = "46213_NVGVSFYADKPEVTQEQK_2", 
#' df_long = example_proteome, example_sample_annotation, 
#' color_by_col = 'm_score', color_by_value = 2)
#'
#' @family feature-level diagnostic functions
#'
#' @export
#' @name feature_level_diagnostics

plot_single_feature  <- function(feature_name, df_long, sample_annotation = NULL,
                                 sample_id_col = 'FullRunName',
                                 measure_col = 'Intensity',
                                 feature_id_col = 'peptide_group_label',
                                 geom = c('point', 'line'),
                                 qual_col = NULL, qual_value = NULL,
                                 batch_col = 'MS_batch',
                                 color_by_batch = FALSE, color_scheme = 'brewer',
                                 order_col = 'order',
                                 vline_color ='red',
                                 facet_col = NULL,
                                 plot_title = NULL,
                                 theme = 'classic',
                                 ylimits = NULL){
  
  #to ensure that missing measurements are NAs (to make them disconnected)
  df_long = df_long %>% complete(!!!syms(c(feature_id_col, sample_id_col)))
  
  #reduce df to measurements of selected features
  plot_df = df_long %>%
    filter(!!(sym(feature_id_col)) %in% feature_name)
  
  #Check the consistency of sample annotation sample IDs and measurement table sample IDs
  plot_df = check_sample_consistency(sample_annotation, sample_id_col, plot_df)
  

  
  #Defining sample order for plotting
  sample_order = define_sample_order(order_col, sample_annotation, facet_col, batch_col, plot_df, 
                                     sample_id_col, color_by_batch)
  order_col = sample_order$order_col
  plot_df = sample_order$df_long
  
  
  
  #Main plotting function
  gg = ggplot(plot_df,
              aes_string(x = order_col, y = measure_col))
  if (identical(geom, 'line')){
    gg = gg + geom_line(color = 'darkgrey', size = .3, aes_string(group = batch_col))
  }
  if (identical(geom, 'point')){
    gg = gg + geom_point()
  }
  if (identical(geom, c('point', 'line'))){
    gg = gg + geom_point() +
      geom_line(color = 'black', alpha = .7, linetype = 'dashed', 
                aes_string(group = batch_col))
  }
  
  #Add coloring for "inferred" measurements / requant values, marked in `color_by_col` with `color_by_value` (e.g. `m_score` and `2`)
  if(!is.null(qual_col)){
    col_data = plot_df %>%
      filter(!!(as.name(qual_col)) == qual_value)
    gg = gg + geom_point(data = col_data,
                         aes_string(x = order_col, y = measure_col),
                         color = 'red', size = 1.5, shape = 8)
  }
  
  #add colors
  gg = color_points_by_batch(color_by_batch, batch_col, gg, color_scheme, sample_annotation)
  if(!is.null(qual_col) && !is.null(color_by_batch) && color_by_batch){
    warning('coloring both inferred values and batches may lead to confusing visualisation, consider plotting separately')
  }
   
  #wrap into facets, if several features are displayed
  #split into facets
  if(!is.null(facet_col)){
    if (length(feature_name) > 1){
      gg = gg + facet_grid(reformulate(facet_col, feature_id_col), 
                           scales = 'free')
    } else {
      gg = gg  + facet_wrap(as.formula(paste("~", facet_col)), 
                            scales = 'free_x')
    }
  } else {
    if (length(feature_name) > 1){
      gg = gg + facet_wrap(as.formula(paste("~", feature_id_col)), 
                           scales = 'free_y')
    }
  }
   
  #add vertical lines, if required (for order-related effects)
  if (!is.null(sample_annotation)){
    gg = add_vertical_batch_borders(order_col, sample_id_col, batch_col, vline_color, 
                                    facet_col, sample_annotation, gg)
  } else {
    gg = add_vertical_batch_borders(order_col, sample_id_col, batch_col, vline_color, 
                                    facet_col, df_long, gg)
  }
  

  
  #Add plot title
  if(!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
    
  #Add the theme
  if (!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic()
  }
  
  #Change the limits of vertical axes
  if(!is.null(ylimits)){
    gg = gg +
      ylim(ylimits)
  }
  
  #Rotate x axis tick labels if the filenames, not numeric order, is displayed
  if (!is.numeric(plot_df[[order_col]])){
    if(is.character(plot_df[[order_col]])){
      plot_df[[order_col]] = factor(plot_df[[order_col]],
                                   levels = unique(plot_df[[order_col]]))
    }
    gg = gg +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  }
  
  #Move the legend to the upper part of the plot to save the horizontal space
  if (length(unique(plot_df[[order_col]])) > 30){
    gg = gg + theme(legend.position="top")
  }
  
  return(gg)
}

#' Plot peptides of one protein
#'
#' Creates a spike-in facetted ggplot2 plot of the value in
#' \code{measure_col} vs \code{order_col} using
#' \code{\link{plot_single_feature}}. Additionally, the resulting plot can also
#' be facetted by batch.
#'
#' @inheritParams feature_level_diagnostics
#' @param protein_name name of the protein as defined in \code{ProteinName}
#' @param ... additional arguments to \code{\link{plot_single_feature}} function
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{spike_ins} containing proteins and (optionally) by \code{batch_col}
#' @examples 
#' peptides_of_one_protein_plot <- plot_peptides_of_one_protein (
#' protein_name = "Haao",  
#' protein_col = "Gene", df_long = example_proteome, 
#' example_sample_annotation, 
#' order_col = 'order', sample_id_col = 'FullRunName', 
#' batch_col = 'MS_batch')
#' 
#' @family feature-level diagnostic functions
#'
#' @export
#'
plot_peptides_of_one_protein <- function(protein_name, peptide_annotation = NULL,
                                         protein_col = 'ProteinName',
                                         df_long, sample_annotation = NULL,
                                         sample_id_col = 'FullRunName',
                                         measure_col = 'Intensity',
                                         feature_id_col = 'peptide_group_label',
                                         geom = c('point', 'line'),
                                         qual_col = NULL, qual_value = NULL,
                                         batch_col = 'MS_batch',
                                         color_by_batch = FALSE, color_scheme = 'brewer',
                                         order_col = 'order',
                                         vline_color ='red',
                                         facet_col = NULL,
                                         plot_title = sprintf('Peptides of %s protein', 
                                                              protein_name),
                                         theme = 'classic', ...){
  
  if (!is.null(peptide_annotation)){
    peptides = peptide_annotation %>%
      filter((!!sym(protein_col)) == protein_name) %>%
      pull(!!sym(feature_id_col)) %>% unique()
    peptides = peptides[peptides %in% df_long[[feature_id_col]]]
  } else {
    peptides = df_long %>%
      filter((!!sym(protein_col)) == protein_name) %>%
      pull(feature_id_col) %>% unique()
  }
  gg = plot_single_feature(peptides, df_long = df_long,
                           sample_annotation = sample_annotation,
                           order_col = order_col,
                           sample_id_col = sample_id_col,
                           batch_col = batch_col, measure_col = measure_col,
                           feature_id_col = feature_id_col,
                           color_by_batch = color_by_batch, 
                           color_scheme = color_scheme,
                           facet_col = facet_col,
                           qual_col = qual_col, 
                           qual_value = qual_value,
                           plot_title = plot_title, ...)
  return(gg)
}

#' Plot spike-in measurements
#'
#' Creates a spike-in facetted ggplot2 plot of the value in
#' \code{measure_col} vs \code{order_col} using
#' \code{\link{plot_single_feature}}. Additionally, the resulting plot can also
#' be facetted by batch.
#'
#' @inheritParams feature_level_diagnostics
#' @param spike_ins substring used to identify spike-in proteins in the column
#'   'ProteinName'
#' @param ... additional arguments to \code{\link{plot_single_feature}} function
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{spike_ins} containing proteins and (optionally) 
#'   by \code{batch_col}
#'
#' @family feature-level diagnostic functions
#' 
#' @examples 
#' spike_in_plot <- plot_spike_in(example_proteome, example_sample_annotation, 
#' protein_col = 'Gene', spike_ins = "BOVINE_A1ag", 
#' plot_title = "Spike-in BOVINE protein peptides")
#' 
#' @export
#'
plot_spike_in <- function(spike_ins = 'BOVIN', peptide_annotation = NULL,
                          protein_col = 'ProteinName',
                          df_long, sample_annotation = NULL,
                          sample_id_col = 'FullRunName',
                          measure_col = 'Intensity',
                          feature_id_col = 'peptide_group_label',
                          geom = c('point', 'line'),
                          qual_col = NULL, qual_value = NULL,
                          batch_col = 'MS_batch',
                          color_by_batch = FALSE, color_scheme = 'brewer',
                          order_col = 'order',
                          vline_color = 'red',
                          facet_col = NULL,
                          plot_title = 'Spike-in BOVINE protein peptides', 
                          theme = 'classic', ...){
  
  if (!is.null(peptide_annotation)){
    df_long = df_long %>%
      merge(peptide_annotation, by = feature_id_col)
  }
  
  if(!is.null(protein_col) && !(protein_col %in% names(df_long))){
    stop('Protein column %s is not found in the data. Check peptide annotation or main data table', protein_col)
  }
  
  if(!is.null(protein_col)){
    spike_in_peptides = df_long %>%
      filter(grepl(spike_ins, !!sym(protein_col))) %>%
      pull(feature_id_col) %>% as.character() %>% unique()
  } else {
    spike_in_peptides = spike_ins
  }
  
  
  gg = plot_single_feature(feature_name = spike_in_peptides, 
                           df_long = df_long,
                           sample_annotation = sample_annotation,
                           sample_id_col = sample_id_col,
                           measure_col = measure_col,
                           feature_id_col = feature_id_col,
                           batch_col = batch_col, 
                           color_by_batch = color_by_batch, 
                           color_scheme = color_scheme,
                           facet_col = facet_col,
                           qual_col = qual_col,
                           qual_value = qual_value,
                           order_col = order_col,
                           plot_title = plot_title, 
                           theme = theme, ...)
  return(gg)
}


#' Plot iRT measurements
#'
#' Creates a iRT facetted ggplot2 plot of the value in
#' \code{measure_col} vs \code{order_col} using
#' \code{\link{plot_single_feature}}. Additionally, the resulting plot can also
#' be facetted by batch.
#'
#' @inheritParams feature_level_diagnostics
#' @param irt_pattern substring used to identify iRT proteins in the column
#'   'ProteinName'
#' @param ... additional arguments to \code{\link{plot_single_feature}} function
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{irt_pattern} containing proteins 
#'   and (optionally) by \code{batch_col}
#'
#' @family feature-level diagnostic functions
#'
#' @examples 
#' irt_plot <- plot_iRT(example_proteome, 
#' example_sample_annotation, 
#' protein_col = 'Gene', irt_pattern = "BOVINE_A1ag")
#' 
#' @export
#'
plot_iRT <- function(irt_pattern = 'iRT',
                     peptide_annotation = NULL,
                     protein_col = 'ProteinName',
                     df_long, sample_annotation = NULL,
                     sample_id_col = 'FullRunName',
                     measure_col = 'Intensity',
                     feature_id_col = 'peptide_group_label',
                     geom = c('point', 'line'),
                     qual_col = NULL, qual_value = NULL,
                     batch_col = 'MS_batch',
                     color_by_batch = FALSE, color_scheme = 'brewer',
                     order_col = 'order',
                     vline_color = 'red',
                     facet_col = NULL,
                     plot_title = 'iRT peptide profile', 
                     theme = 'classic', ...){
  
  if (!is.null(peptide_annotation)){
    df_long = df_long %>%
      merge(peptide_annotation, by = feature_id_col)
  }
  iRT_peptides = df_long %>%
    filter(grepl(irt_pattern, !!sym(protein_col))) %>%
    pull(feature_id_col)  %>% unique()
  gg = plot_single_feature(feature_name = iRT_peptides, 
                           df_long = df_long, 
                           sample_annotation = sample_annotation,
                           sample_id_col = sample_id_col,
                           measure_col = measure_col,
                           feature_id_col = feature_id_col,
                           geom = geom,
                           qual_col = qual_col, qual_value = qual_value,
                           batch_col = batch_col,
                           color_by_batch = color_by_batch, color_scheme = color_scheme,
                           order_col = order_col, 
                           vline_color = vline_color,
                           facet_col = facet_col,
                           plot_title = plot_title, 
                           theme = theme, ...)
  return(gg)
}


#' Plot peptide measurements across multi-step analysis
#'
#' Plot Intensity of a few representative peptides for each step of the analysis
#' including the fitting curve
#'
#' @inheritParams feature_level_diagnostics
#' @param fit_df data frame typically output generated from nonlinear curve 
#'   fitting by \code{normalize_custom_fit}
#' @param fit_value_var column denoting intensity values, typically fitted to curve
#' @param ... additional arguments to \code{\link{plot_single_feature}} function
#'
#' @return \code{ggplot}-class plot with minimally two facets (before and after
#'   non-linear fit) with \code{measure_col} (Intensity) vs \code{order_col}
#'   (injection order) for selected peptides (specified in \code{feature_name})
#'
#' @family feature-level diagnostic functions
#' @examples 
#' loess_fit_70 <- adjust_batch_trend(example_proteome, 
#' example_sample_annotation, span = 0.7)
#' 
#' fitting_curve_plot <- plot_with_fitting_curve(feature_name = "10231_QDVDVWLWQQEGSSK_2", 
#' df_long = example_proteome, example_sample_annotation, 
#' fit_df = loess_fit_70$fit_df, plot_title = "Curve fitting with 70% span")
#'
#' @export
#'

plot_with_fitting_curve <- function(feature_name, 
                                    fit_df, fit_value_var = 'fit', 
                                    df_long, sample_annotation = NULL,
                                    sample_id_col = 'FullRunName',
                                    measure_col = 'Intensity',
                                    feature_id_col = 'peptide_group_label',
                                    geom = c('point', 'line'),
                                    qual_col = NULL, qual_value = NULL, 
                                    batch_col = 'MS_batch',
                                    color_by_batch = FALSE, color_scheme = 'brewer', 
                                    order_col = 'order',
                                    vline_color = 'grey',
                                    facet_col = NULL,
                                    plot_title = sprintf("Fitting curve of %s peptide", 
                                                         feature_name),
                                    theme = 'classic'){
  
  if(length(feature_name) > 10){
    warning("Visualisation of individual features can be suboptimal,
            consider exploring no more than 5 features at a time")
  }
  #Plotting single features as usually (only batch coloring, if specified, is on fitting-curve layer)
  gg = plot_single_feature(feature_name = feature_name, df_long = df_long,
                           sample_annotation = sample_annotation,
                           sample_id_col = sample_id_col,
                           measure_col = measure_col,
                           feature_id_col = feature_id_col,
                           geom = geom,
                           qual_col = qual_col, qual_value =  qual_value,
                           batch_col = batch_col,
                           color_by_batch = FALSE, color_scheme = NULL,
                           order_col = order_col,
                           vline_color = vline_color, 
                           facet_col = facet_col,
                           plot_title = plot_title, 
                           theme = theme)
  
  fit_df = fit_df %>%
    filter(!!(sym(feature_id_col)) %in% feature_name)
  
  if (!is.null(sample_annotation)){
    fit_df = fit_df %>%
      merge(sample_annotation, by = c(sample_id_col))
  }
   
  if(identical(color_by_batch, FALSE)){
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_var, x = order_col, 
                                   group = batch_col), 
                        color = 'red')
  } else {
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_var, x = order_col,
                                   group = batch_col, 
                                   color = batch_col), size = 1.25)
    #TODO: add coloring scheme as in other functions here, too
    if(color_by_batch & length(color_scheme) == length(unique(fit_df[[batch_col]]))){
      gg = gg + scale_color_manual(values = color_scheme)
    }
  }
  
  return(gg)
}
