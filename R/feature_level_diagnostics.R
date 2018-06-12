#' Plot peptide measurements
#'
#' Creates a peptide facetted ggplot2 plot of the value in \code{measure_col}
#' vs \code{order_col}. Additionally, the resulting plot can also be facetted
#' by batch.
#'
#' @inheritParams proBatch
#' @param pep_name name of the peptide for diagnostic profiling
#' @param geom whether to show the feature as points and/or connect by lines
#' @param color_by_batch (logical) whether to color points by batch
#' @param facet_by_batch (logical) whether to plot each batch in its own facet
#' @param plot_title the string indicating the source of the peptides
#' @param requant if data frame: requant values; if logical: whether to indicate
#'   requant values (requires 'requant' or 'm_score' column in \code{df_long})
#' @param theme plot theme (default is 'classical'; other options not
#'   implemented)
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{pep_name} and (optionally) by \code{batch_col}
#'
#' @family feature-level diagnostic functions
#'
#' @export
#'
plot_peptide_level <- function(pep_name, df_long, sample_annotation,
                               order_col = NULL,
                               sample_id_col = 'FullRunName',
                               batch_col = 'MS_batch.final',
                               measure_col = 'Intensity',
                               feature_id_col = 'peptide_group_label',
                               geom = c('point', 'line'),
                               color_by_batch = F, facet_by_batch = F,
                               plot_title = NULL, requant = NULL,
                               theme = 'classic'){
  #TODO: suggest faceting by instrument
  #TODO: plot fit after LOESS



  plot_df = df_long %>%
    filter(UQ(as.name(feature_id_col)) %in% pep_name)
  if (!all(names(sample_annotation) %in% names(df_long))){
    sample_annotation = sample_annotation %>% arrange_(order_col)
    #remove annotation columns, except the sample_id_col
    common_cols = intersect(names(sample_annotation), names(plot_df))
    cols_to_remove = setdiff(common_cols, sample_id_col)
    plot_df = plot_df %>%
      select(-one_of(cols_to_remove))
    plot_df = plot_df %>%
      merge(sample_annotation, by = sample_id_col)
  }


  if(is.null(order_col)){
    warning("order column wasn't specified, putting row number as an order within a batch")
    plot_df = plot_df %>%
      group_by_at(vars(one_of(batch_col))) %>%
      mutate(order = row_number())
    order_col = 'order'
  }



  gg = ggplot(plot_df,
              aes_string(x = order_col, y = measure_col))
  if (identical(geom, 'line')){
    gg = gg + geom_line(color = 'darkgrey', size = .3)
  }
  if (identical(geom, 'point')){
    gg = gg + geom_point()
  }

  if (identical(geom, c('point', 'line'))){
    gg = gg + geom_point() +
      geom_line(color = 'black', alpha = .7, linetype = 'dashed')
  }
  if(!color_by_batch & !is.null(batch_col)){
    batch.tipping.points = cumsum(table(sample_annotation[[batch_col]]))+.5
    gg = gg + geom_vline(xintercept = batch.tipping.points,
                         color = 'grey', linetype = 'dashed')
  } else {
    gg = gg + geom_point(aes_string(color = batch_col))
    if(facet_by_batch){
      if (length(pep_name) > 1){
        gg = gg + facet_grid(reformulate(batch_col, pep_name), scales = 'free_y')
      } else {
        gg = gg  + facet_wrap(as.formula(paste("~", batch_col)), scales = 'free_y')
      }
    }
  }
  if (length(pep_name) > 1){
    gg = gg + facet_wrap(as.formula(paste("~", feature_id_col)), scales = 'free_y')
  }
  if(!is.null(plot_title)){
    gg = gg + ggtitle(plot_title)
  }
  if (!is.null(requant)){
    if (is.data.frame(requant)){
      data_requant = requant
    } else {
      if('requant' %in% names(df_long)){
        data_requant = plot_df %>% filter(requant)
      } else {
        if('m_score' %in% names(df_long)){
          data_requant = plot_df %>% filter(m_score == 2)
        } else{
          stop('requant cannot be plotted without requant table or m_scores!')
        }
      }
    }
    data_requant = data_requant %>%
      filter(UQ(as.name(feature_id_col)) %in% pep_name)
    gg = gg + geom_point(data = data_requant,
                         aes_string(x = order_col, y = measure_col),
                         color = 'red', size = .3, shape = 8)
  }
  if (theme == 'classic'){
    gg = gg + theme_classic()
  }
  return(gg)
}

#' Plot spike-in measurements
#'
#' Creates a spike-in facetted ggplot2 plot of the value in
#' \code{measure_col} vs \code{order_col} using
#' \code{\link{plot_peptide_level}}. Additionally, the resulting plot can also
#' be facetted by batch.
#'
#' @inheritParams plot_peptide_level
#' @param spike_ins substring used to identify spike-in proteins in the column
#'   'ProteinName'
#' @param ... additional arguments to \code{\link{plot_peptide_level}} function
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{spike_ins} containing proteins and (optionally) by \code{batch_col}
#'
#' @family feature-level diagnostic functions
#'
#' @export
#'
plot_spike_ins <- function(df_long, sample_annotation,
                           order_col = 'order',
                           spike_ins = 'BOVIN',
                           sample_id_col = 'FullRunName',
                           batch_col = 'MS_batch',
                           measure_col = 'Intensity',
                           feature_id_col = 'peptide_group_label',
                           plot_title = 'Spike-in BOVINE protein peptides', ...){
  spike_in_peptides = df_long %>%
    filter(grepl(spike_ins, ProteinName)) %>%
    pull(feature_id_col) %>% unique()
  gg = plot_peptide_level(spike_in_peptides, df_long = df_long,
                          sample_annotation = sample_annotation,
                          order_col = order_col,
                          sample_id_col = sample_id_col,
                          batch_col = batch_col, measure_col = measure_col,
                          feature_id_col = feature_id_col,
                          plot_title = plot_title, ...)
  return(gg)
}


#' Plot iRT measurements
#'
#' Creates a iRT facetted ggplot2 plot of the value in
#' \code{measure_col} vs \code{order_col} using
#' \code{\link{plot_peptide_level}}. Additionally, the resulting plot can also
#' be facetted by batch.
#'
#' @inheritParams plot_peptide_level
#' @param irt_pattern substring used to identify irts proteins in the column
#'   'ProteinName'
#' @param ... additional arguments to \code{\link{plot_peptide_level}} function
#'
#' @return ggplot2 type plot of \code{measure_col} vs \code{order_col},
#'   faceted by \code{irt_pattern} containing proteins and (optionally) by \code{batch_col}
#'
#' @family feature-level diagnostic functions
#'
#' @export
#'
#' @examples
plot_iRTs <- function(df_long, sample_annotation,
                      order_col = NULL,
                      irt_pattern = 'iRT',
                      batch_col = 'MS_batch.final',
                      sample_id_col = 'FullRunName',
                      feature_id_col = 'peptide_group_label',
                      measure_col = 'Intensity',
                      plot_title = 'iRT peptide profile', ...){
  iRT_peptides = df_long %>%
    filter(grepl(irt_pattern, ProteinName)) %>%
    pull(feature_id_col)  %>% unique()
  gg = plot_peptide_level(iRT_peptides, df_long, sample_annotation,
                          order_col = order_col,
                          sample_id_col = sample_id_col,
                          batch_col = batch_col,
                          feature_id_col = feature_id_col,
                          measure_col = measure_col,
                          plot_title = plot_title, ...)
  return(gg)
}


#' Plot peptide measurements across multi-step analysis
#'
#' Plot Intensity of a few representative peptides for each step of the analysis
#' including the fitting curve
#'
#' @inheritParams plot_peptide_level
#' @param pep_name name of the peptide for diagnostic profiling
#' @param data_df_all_steps data frame, similar to \code{df_long}
#'   \link{proBatch},  where each row is a single feature in a single sample, at
#'   a certain step of the analysis (minimally raw and after linear
#'   normalization) thus it has minimally the following columns:
#'   \code{sample_id_col}, \code{feature_id_col}, \code{measure_col}, and
#'   \code{fit_step}, but usually also \code{m_score}
#' @param fit_df
#' @param fit_value_var
#' @param geom for the intensity \code{measure_col} profile:
#'
#' @return \code{ggplot}-class plot with minimally two facets (before and after
#'   non-linear fit) with \code{measure_col} (Intensity) vs \code{order_col}
#'   (injection order) for selected peptides (specified in \code{pep_name})
#'
#' @family feature-level diagnostic functions
#'
#' @export
#'

# TODO: Add descriptions of fit_df and fit_value_var
plot_with_fitting_curve <- function(pep_name, data_df_all_steps,
                                    sample_annotation,
                                    fit_df,
                                    fit_value_var = 'fit', fit_step = '3_loess_fit',
                                    order_col = NULL,
                                    sample_id_col = 'FullRunName',
                                    batch_col = 'MS_batch',
                                    measure_col = 'Intensity',
                                    feature_id_col = 'peptide_group_label',
                                    geom = c('point', 'line'),
                                    color_by_batch = F, facet_by_batch = F,
                                    plot_title = NULL, requant = NULL,
                                    theme = 'classic'){
  #TODO: in non-linear fit, change "Intensity_normalized" to "Intensity", but rename "Intensity" as "Intensity_before_fit"
  if(length(pep_name) > 10){
    warning("Visualisation of individual features can be suboptimal,
            consider exploring no more than 5 features at a time")
  }
  gg = plot_peptide_level(pep_name, df_long = data_df_all_steps,
                          sample_annotation = sample_annotation,
                          order_col = order_col,
                          sample_id_col = sample_id_col,
                          batch_col = batch_col,
                          measure_col = measure_col,
                          feature_id_col = feature_id_col,
                          plot_title = plot_title,
                          facet_by_batch = facet_by_batch)
  if(!("Step" %in% names(fit_df))){
    fit_df$Step = fit_step
    if (!(fit_step %in% data_df_all_steps$Step)){
      stop('specify step for which the curve fit should be shown')
    }
  }
  fit_df = fit_df %>%
    filter(UQ(as.name(feature_id_col)) %in% pep_name) %>%
    merge(sample_annotation, by = c(sample_id_col, batch_col))
  if(identical(color_by_batch, FALSE)){
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_var, x = order_col, group = batch_col),
                        color = 'red')+
      facet_grid(as.formula(paste(c(feature_id_col, "~", 'Step'), collapse = ' ')),
                 scales = 'free_y')
  } else {
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_var, x = order_col,
                                   group = batch_col, color = batch_col), size = 1.25)+
      facet_grid(as.formula(paste(c(feature_id_col, "~", 'Step'), collapse = ' ')),
                 scales = 'free_y')
    if(length(color_by_batch) == length(unique(fit_df[[batch_col]]))){
      gg = gg + scale_color_manual(values = color_by_batch)
    }
  }

  return(gg)
}
