#' @title Ploting peptide measurements
#'
#' @description Creates a peptide faceted ggplot2 plot of the value in 
#' \code{measure_col}
#' vs \code{order_col} (if `NULL`, x-axis is simply a sample name order). 
#' Additionally, the resulting plot can also be colored either by batch factor, 
#' by quality factor (e.g. imputated/non-imputed) and, if needed, faceted by 
#' another batch factor, e.g. an instrument.
#'  If the non-linear curve was fit, this can also be added to the plot, see 
#'  functions specific to each case below
#'
#' @inheritParams proBatch
#' @param feature_name name of the selected feature (e.g. peptide) for 
#' diagnostic profiling
#' @param geom whether to show the feature as points and/or connect by lines 
#' (accepted values are: 1. \code{point}, \code{line} and 
#' \code{c('point', 'line')})
#' @param protein_name name of the protein as defined in \code{ProteinName}
#' @param irt_pattern substring used to identify iRT proteins in the column
#'   'ProteinName'
#' @param spike_ins name of feature(s), typically proteins that were spiked in 
#' for control
#' @param vline_color color of vertical lines, typically separating 
#'  different MS batches in ordered runs; 
#'  should be `NULL` for experiments without intrinsic order
#' @param ylimits range of y-axis to plot feature-level trends 
#' @param fit_df data frame output of \code{adjust_batch_trend_df} to be plotted
#' with the line
#' @param fit_value_col column in \code{fit_df} where the values for fitting 
#' trend are found
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{feature_name} and (optionally) by \code{batch_col}
#' @examples 
#' single_feature_plot <- plot_single_feature(feature_name = c("46213_NVGVSFYADKPEVTQEQK_2","10081_NVQGIIDILK_2"), 
#' df_long = example_proteome, example_sample_annotation, 
#' qual_col = NULL)
#' 
#' #color measurements by factor, related to order (MS_batch)
#' plot_single_feature(feature_name = c("46213_NVGVSFYADKPEVTQEQK_2","10081_NVQGIIDILK_2"), 
#' df_long = example_proteome, example_sample_annotation, 
#' qual_col = NULL, color_by_batch = TRUE, batch_col = 'MS_batch')
#' 
#' #color measurements by factor, with order-unrelated factor
#' single_feature_plot <- plot_single_feature(feature_name = c("46213_NVGVSFYADKPEVTQEQK_2","10081_NVQGIIDILK_2"), 
#' df_long = example_proteome, example_sample_annotation, 
#' qual_col = NULL, color_by_batch = TRUE, batch_col = 'Diet', geom = 'point', 
#' vline_color = NULL)
#' 
#' #saving the plot
#' \dontrun{
#' single_feature_plot <- plot_single_feature(feature_name = c("46213_NVGVSFYADKPEVTQEQK_2","10081_NVQGIIDILK_2"), 
#' df_long = example_proteome, example_sample_annotation, 
#' qual_col = NULL, filename = 'test_peptide.png', 
#' width = 28, height = 18, units = 'cm')
#' }
#' 
#' #to examine peptides of a single protein:
#' peptides_of_one_protein_plot <- plot_peptides_of_one_protein (
#' protein_name = "Haao", peptide_annotation = example_peptide_annotation,
#' protein_col = "Gene", df_long = example_proteome, 
#' sample_annotation = example_sample_annotation, 
#' order_col = 'order', sample_id_col = 'FullRunName', 
#' batch_col = 'MS_batch')
#' 
#' #saving the peptides of one protein
#' \dontrun{
#'  peptides_of_one_protein_plot <- plot_peptides_of_one_protein (
#' protein_name = "Haao", peptide_annotation = example_peptide_annotation,
#' protein_col = "Gene", df_long = example_proteome, 
#' sample_annotation = example_sample_annotation, 
#' order_col = 'order', sample_id_col = 'FullRunName', 
#' batch_col = 'MS_batch',
#' filename = 'test_protein.png', width = 14, height = 9, units = 'in')}
#' 
#' #to illustrate spike-ins:
#' spike_in_plot <- plot_spike_in(spike_ins = "BOVINE_A1ag", 
#' peptide_annotation = example_peptide_annotation, protein_col = 'Gene', 
#' df_long = example_proteome, sample_annotation = example_sample_annotation, 
#' sample_id_col = 'FullRunName',
#' plot_title = "Spike-in BOVINE protein peptides")
#' 
#' #to illustrate iRT peptides:
#' irt_plot <- plot_iRT(irt_pattern = "iRT", 
#' peptide_annotation = example_peptide_annotation, 
#' df_long = example_proteome, sample_annotation = example_sample_annotation, 
#' protein_col = 'Gene')
#'
#' #illustrate the fitting curve:
#' special_peptide = example_proteome$peptide_group_label == "10231_QDVDVWLWQQEGSSK_2"
#' loess_fit_70 <- adjust_batch_trend_df(example_proteome[special_peptide,], 
#' example_sample_annotation, span = 0.7)
#' 
#' fitting_curve_plot <- plot_with_fitting_curve(feature_name = "10231_QDVDVWLWQQEGSSK_2", 
#' df_long = example_proteome, sample_annotation = example_sample_annotation, 
#' fit_df = loess_fit_70, plot_title = "Curve fitting with 70% span")
#' 
#' #with curves colored by the corresponding batch:
#' fitting_curve_plot <- plot_with_fitting_curve(feature_name = "10231_QDVDVWLWQQEGSSK_2", 
#' df_long = example_proteome, sample_annotation = example_sample_annotation, 
#' fit_df = loess_fit_70, plot_title = "Curve fitting with 70% span", 
#' color_by_batch = TRUE, batch_col = 'MS_batch')
#'
#' @name feature_level_diagnostics
NULL

#'
#' @export
#' @rdname feature_level_diagnostics
plot_single_feature  <- function(feature_name, df_long, 
                                 sample_annotation = NULL,
                                 sample_id_col = 'FullRunName',
                                 measure_col = 'Intensity',
                                 feature_id_col = 'peptide_group_label',
                                 geom = c('point', 'line'),
                                 qual_col = NULL, qual_value = NULL,
                                 batch_col = 'MS_batch',
                                 color_by_batch = FALSE, 
                                 color_scheme = 'brewer',
                                 order_col = 'order',
                                 vline_color ='red',
                                 facet_col = NULL,
                                 filename = NULL, width = NA, height = NA, 
                                 units = c('cm','in','mm'),
                                 plot_title = NULL,
                                 theme = 'classic',
                                 ylimits = NULL){
  
  #to ensure that missing measurements are NAs (to make them disconnected)
  #df_long = df_long %>% complete(!!!syms(c(feature_id_col, sample_id_col)))
  
  #reduce df to measurements of selected features
  plot_df = df_long %>%
    filter(!!(sym(feature_id_col)) %in% feature_name)
  rm(df_long)
  gc()
  
  #Check the consistency of sample annot. sample IDs and measur.table sample IDs
  plot_df = check_sample_consistency(sample_annotation, sample_id_col, plot_df, 
                                     batch_col, order_col, facet_col)
  

  
  #Defining sample order for plotting
  sample_order = define_sample_order(order_col, sample_annotation, facet_col, 
                                     batch_col, plot_df, 
                                     sample_id_col, color_by_batch)
  order_col = sample_order$order_col
  plot_df = sample_order$df_long
  
  #Ensure that batch-coloring-related arguments are defined properly
  if(!is.null(batch_col)){
    if(!(batch_col %in% names(plot_df))){
      stop('batches cannot be colored as the batch column or sample ID column
           is not defined, check sample_annotation and data matrix')
    }
    } else {
      if (color_by_batch){
        warning('batches cannot be colored as the batch column is defined as 
                NULL, continuing without colors')
        color_by_batch = FALSE
      }
    }
  
  #For order definition and subsequent faceting, facet column has to be in the 
  #data frame
  if(!is.null(facet_col)){
    if ( !(facet_col %in% names(plot_df))){
      stop(sprintf('"%s" is specified as column for faceting, but is not present 
                    in the data, check sample annotation data frame', 
                   facet_col))
    }
  }
  
  if (!is.null(batch_col)){
    batch_vector <- sample_annotation[[batch_col]]
    is_factor = is_batch_factor(batch_vector, color_scheme)
  }
  
  #Main plotting function
  gg = ggplot(plot_df,
              aes_string(x = order_col, y = measure_col))
  if (identical(geom, 'line')){
    gg = gg + geom_line(color = 'darkgrey', size = .3, 
                        aes_string(group = batch_col))
  }
  if (identical(geom, 'point')){
    gg = gg + geom_point()
  }
  if (identical(geom, c('point', 'line'))){
    if (is.null(batch_col) || !is_factor){
      gg = gg + geom_point() +
        geom_line(color = 'black', alpha = .7, linetype = 'dashed')
    } else {
      if(is_factor){
        gg = gg + geom_point() +
          geom_line(color = 'black', alpha = .7, linetype = 'dashed', 
                    aes_string(group = batch_col))
      } 
    }
  }
  
  #Add coloring for "inferred" measurements / imputed (requant) values, marked 
  #in `color_by_col` with `color_by_value` (e.g. `m_score` and `2`)
  if(!is.null(qual_col)){
    col_data = plot_df %>%
      filter(!!(as.name(qual_col)) == qual_value)
    gg = gg + geom_point(data = col_data,
                         aes_string(x = order_col, y = measure_col),
                         color = 'red', size = 1, shape = 8)
  }
  
  
  
  #add colors
  gg = color_by_factor(color_by_batch = color_by_batch, 
                       batch_col = batch_col, gg = gg, 
                       color_scheme = color_scheme, 
                       sample_annotation = plot_df,
                       fill_or_color = 'color')
  if(!is.null(qual_col) && !is.null(color_by_batch) && color_by_batch){
    warning('coloring both inferred values and batches may lead to confusing 
            visualisation, consider plotting separately')
  }
   
  #wrap into facets, if several features are displayed
  #split into facets
  if(!is.null(facet_col)){
    if (facet_col != feature_id_col && length(feature_name) > 1){
      gg = gg + facet_grid(reformulate(facet_col, feature_id_col), 
                           scales = 'free')
    } else {
      if (facet_col == feature_id_col && length(feature_name) >1){
        gg = gg + facet_wrap(as.formula(paste("~", feature_id_col)), 
                             scales = 'free_y')
      } else {
        gg = gg  + facet_wrap(as.formula(paste("~", facet_col)), 
                              scales = 'free_x')
      }
    }
  } else {
    if (length(feature_name) > 1){
      gg = gg + facet_wrap(as.formula(paste("~", feature_id_col)), 
                           scales = 'free_y')
    }
  }
  
  #add vertical lines, if required (for order-related effects)
  if (!is.null(batch_col) && is_factor){
    gg = add_vertical_batch_borders(order_col, sample_id_col, batch_col, 
                                    vline_color, 
                                    facet_col, plot_df, gg)
  }
  
  #Add plot title
  if(!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
    
  #Add the theme
  if (!is.null(theme) && theme == 'classic'){
    gg = gg + theme_classic()
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
  if (!is.numeric(plot_df[[order_col]])){
    if(is.character(plot_df[[order_col]])){
      plot_df[[order_col]] = factor(plot_df[[order_col]],
                                   levels = unique(plot_df[[order_col]]))
    }
    gg = gg +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
  }
  
  #Move the legend to the upper part of the plot to save the horizontal space
  if (length(unique(plot_df[[order_col]])) > 30 && color_by_batch && is_factor){
    gg = gg + theme(legend.position="top")
  }
  
  #save the plot
  save_ggplot(filename, units, width, height, gg)
  
  return(gg)
}

#'
#' @export
#' @rdname feature_level_diagnostics
#'
plot_peptides_of_one_protein <- function(protein_name, 
                                         peptide_annotation = NULL,
                                         protein_col = 'ProteinName',
                                         df_long, sample_annotation = NULL,
                                         sample_id_col = 'FullRunName',
                                         measure_col = 'Intensity',
                                         feature_id_col = 'peptide_group_label',
                                         geom = c('point', 'line'),
                                         qual_col = NULL, qual_value = NULL,
                                         batch_col = 'MS_batch',
                                         color_by_batch = FALSE, 
                                         color_scheme = 'brewer',
                                         order_col = 'order',
                                         vline_color ='red',
                                         facet_col = NULL,
                                         filename = NULL, 
                                         width = NA, height = NA, 
                                         units = c('cm','in','mm'),
                                         plot_title = sprintf('Peptides of %s protein', 
                                                              protein_name),
                                         theme = 'classic'){
  
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
                           sample_id_col = sample_id_col,
                           measure_col = measure_col,
                           feature_id_col = feature_id_col,
                           geom = geom,
                           qual_col = qual_col, 
                           qual_value = qual_value,
                           batch_col = batch_col, 
                           color_by_batch = color_by_batch, 
                           color_scheme = color_scheme,
                           order_col = order_col,
                           vline_color = vline_color,
                           facet_col = facet_col,
                           plot_title = plot_title, 
                           theme = theme)
  
  #save the plot
  save_ggplot(filename, units, width, height, gg)
  return(gg)
}


#' 
#' @export
#' @rdname feature_level_diagnostics
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
                          filename = NULL, width = NA, height = NA, 
                          units = c('cm','in','mm'),
                          plot_title = sprintf('Spike-in %s plots', spike_ins), 
                          theme = 'classic'){
  
  if(!is.null(protein_col)){
    if(protein_col %in% names(df_long)){
      spike_in_peptides = df_long %>%
        filter(grepl(spike_ins, !!sym(protein_col))) %>%
        pull(feature_id_col) %>% as.character() %>% unique()
    } else{
      if (!is.null(peptide_annotation)){
        if(protein_col %in% names(peptide_annotation)){
          spike_in_peptides = peptide_annotation %>%
            filter(grepl(spike_ins, !!sym(protein_col))) %>%
            pull(feature_id_col) %>% as.character() %>% unique()
          df_long = df_long %>%
            filter(!!sym(feature_id_col) %in% spike_in_peptides) %>%
            inner_join(peptide_annotation, by = feature_id_col)
        }
        
      }
    }
  } else {
    spike_in_peptides = spike_ins
  }
  
  
  
  if(!is.null(protein_col) && !(protein_col %in% names(df_long))){
    stop('Protein column %s is not found in the data. Check peptide annotation 
         or main data table', protein_col)
  }
  
  
  
  
  gg = plot_single_feature(feature_name = spike_in_peptides, 
                           df_long = df_long,
                           sample_annotation = sample_annotation,
                           sample_id_col = sample_id_col,
                           measure_col = measure_col,
                           feature_id_col = feature_id_col,
                           geom = geom,
                           qual_col = qual_col, qual_value = qual_value,
                           batch_col = batch_col, 
                           color_by_batch = color_by_batch, 
                           color_scheme = color_scheme,
                           order_col = order_col,
                           facet_col = facet_col,
                           plot_title = plot_title, 
                           theme = theme)
  
  #save the plot
  save_ggplot(filename, units, width, height, gg)
  return(gg)
}

#' 
#' @export
#' @rdname feature_level_diagnostics
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
                     filename = NULL, width = NA, height = NA, 
                     units = c('cm','in','mm'),
                     plot_title = 'iRT peptide profile', 
                     theme = 'classic'){
  
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
                           color_by_batch = color_by_batch, 
                           color_scheme = color_scheme,
                           order_col = order_col, 
                           vline_color = vline_color,
                           facet_col = facet_col,
                           plot_title = plot_title, 
                           theme = theme)
  
  #save the plot
  save_ggplot(filename, units, width, height, gg)
  return(gg)
}

#'
#' @export
#' @rdname feature_level_diagnostics

plot_with_fitting_curve <- function(feature_name, 
                                    fit_df, fit_value_col = 'fit', 
                                    df_long, sample_annotation = NULL,
                                    sample_id_col = 'FullRunName',
                                    measure_col = 'Intensity',
                                    feature_id_col = 'peptide_group_label',
                                    geom = c('point', 'line'),
                                    qual_col = NULL, qual_value = NULL, 
                                    batch_col = 'MS_batch',
                                    color_by_batch = FALSE, 
                                    color_scheme = 'brewer', 
                                    order_col = 'order',
                                    vline_color = 'grey',
                                    facet_col = NULL,
                                    filename = NULL, width = NA, height = NA, 
                                    units = c('cm','in','mm'),
                                    plot_title = sprintf("Fitting curve of %s 
                                                         peptide", 
                                                         paste(feature_name, 
                                                               collapse = ' ')),
                                    theme = 'classic'){
  
  if(length(feature_name) > 10){
    warning("Visualisation of individual features can be suboptimal,
            consider exploring no more than 5 features at a time")
  }
  #Plotting single features as usually (only batch coloring, if specified, is on
  #fitting-curve layer)
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
    fit_df = check_sample_consistency(sample_annotation = sample_annotation, 
                                      df_long = fit_df,
                                      sample_id_col = sample_id_col, 
                                      batch_col = batch_col, 
                                      order_col = order_col, facet_col = facet_col)
    fit_df = define_sample_order(order_col = order_col, 
                                 sample_annotation = sample_annotation,
                                 df_long = fit_df, 
                                 sample_id_col = sample_id_col, 
                                 facet_col = facet_col, batch_col = batch_col, 
                                 color_by_batch = color_by_batch)$df_long
  }
   
  if(color_by_batch && !is.null(batch_col)){
    
    batch_vector <- sample_annotation[[batch_col]]
    n_batches <- length(unique(batch_vector))
    is_factor = is_batch_factor(batch_vector, color_scheme)
    if(!is_factor){
      stop('coloring by fitting curve possible only for the batch factors 
           corresponding to curve-fitting batches. Change color_by_batch = F 
           to see the curves without colors or batch_col to the right 
           batch factor')
    }
    
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_col, x = order_col,
                                   group = batch_col, 
                                   color = batch_col), size = 1.25)
    
    gg = add_color_scheme_discrete(color_scheme, n_batches, 
                                   fill_or_color = 'color', 
                                   gg = gg, batch_col = batch_col)
  } else {
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_col, x = order_col, 
                                   group = batch_col), 
                        color = 'red', size = 1.25)
  }
  
  #save the plot
  save_ggplot(filename, units, width, height, gg)
  
  return(gg)
}
