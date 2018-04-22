#' plot a single peptide or several peptides each in its own facet
#'
#' @param pep_name
#' @param data_df_long
#' @param sample_annotation
#' @param order_column
#' @param measurement.col
#' @param batch_column
#' @param feature_id_column
#' @param geom
#' @param facet_by_batch
#' @param title
#' @param requant
#' @param sample_id_column
#' @param color_by_batch
#' @param theme
#'
#' @return
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @examples
plot_peptide_level <- function(pep_name, data_df_long, sample_annotation,
                               batch_column = 'MS_batch.final',
                               feature_id_column = 'peptide_group_label',
                               measurement.col = 'Intensity',
                               sample_id_column = 'FullRunName',
                               order_column = NULL, geom = c('point', 'line'),
                               color_by_batch = F, facet_by_batch = F,
                               title = NULL, requant = NULL, theme = 'classic'){
  #TODO: suggest faceting by instrument
  #TODO: plot fit after LOESS



  plot_df = data_df_long %>%
    filter(rlang::UQ(as.name(feature_id_column)) %in% pep_name)
  if (!all(names(sample_annotation) %in% names(data_df_long))){
    sample_annotation = sample_annotation %>% arrange_(order_column)
    #remove annotation columns, except the sample_id_column
    common_cols = intersect(names(sample_annotation), names(plot_df))
    cols_to_remove = setdiff(common_cols, sample_id_column)
    plot_df = plot_df %>%
      select(-one_of(cols_to_remove))
    plot_df = plot_df %>%
      merge(sample_annotation, by = sample_id_column)
  }


  if(is.null(order_column)){
    warning("order column wasn't specified, putting row number as an order within a batch")
    plot_df = plot_df %>%
      group_by_at(vars(one_of(batch_column))) %>%
      mutate(order = row_number())
    order_column = 'order'
  }



  gg = ggplot(plot_df,
              aes_string(x = order_column, y = measurement.col))
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
  if(!color_by_batch & !is.null(batch_column)){
    batch.tipping.points = cumsum(table(sample_annotation[[batch_column]]))+.5
    gg = gg + geom_vline(xintercept = batch.tipping.points,
                         color = 'grey', linetype = 'dashed')
  } else {
    gg = gg + geom_point(aes_string(color = batch_column))
    if(facet_by_batch){
      if (length(pep_name) > 1){
        gg = gg + facet_grid(reformulate(batch_column, pep_name), scales = 'free_y')
      } else {
        gg = gg  + facet_wrap(as.formula(paste("~", batch_column)), scales = 'free_y')
      }
    }
  }
  if (length(pep_name) > 1){
    gg = gg + facet_wrap(as.formula(paste("~", feature_id_column)), scales = 'free_y')
  }
  if(!is.null(title)){
    gg = gg + ggtitle(title)
  }
  if (!is.null(requant)){
    if (is.data.frame(requant)){
      data_requant = requant
    } else {
      if('requant' %in% names(data_df_long)){
        data_requant = plot_df %>% filter(requant)
      } else {
        if('m_score' %in% names(data_df_long)){
          data_requant = plot_df %>% filter(m_score == 2)
        } else{
          stop('requant cannot be plotted without requant table or m_scores!')
        }
      }
    }
    data_requant = data_requant %>%
      filter(rlang::UQ(as.name(feature_id_column)) %in% pep_name)
    gg = gg + geom_point(data = data_requant,
                         aes_string(x = order_column, y = measurement.col),
                         color = 'red', size = .3, shape = 8)
  }
  if (theme == 'classic'){
    gg = gg + theme_classic()
  }
  return(gg)
}

#' Plot iRT peptides
#'
#' @param data_df_long - "openSWATH" format data frame
#' @param sample_annotation
#' @param batch_column
#' @param feature_id_column
#' @param measurement.col
#' @param order_column
#' @param sample_id_column
#' @param ... additional arguments to plot_peptide_level function
#'
#' @return
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
plot_iRTs <- function(data_df_long, sample_annotation,
                      batch_column = 'MS_batch.final',
                      sample_id_column = 'FullRunName',
                      feature_id_column = 'peptide_group_label',
                      measurement.col = 'Intensity',
                      order_column = 'order',  ...){
  iRT_peptides = data_df_long %>%
    filter(grepl('iRT', ProteinName)) %>%
    pull(feature_id_column)  %>% unique()
  gg = plot_peptide_level(iRT_peptides, data_df_long, sample_annotation,
                          batch_column = batch_column,
                          feature_id_column = feature_id_column,
                          measurement.col = measurement.col,
                          sample_id_column = sample_id_column,
                          order_column = order_column,
                          title = 'iRT_peptides', ...)
  return(gg)
}

#' Plot Spike-in peptides/proteins
#'
#' @param data_df_long
#' @param sample_annotation
#' @param spike_ins
#' @param order_column
#' @param measurement.col
#' @param batch_column
#' @param sample_id_column
#' @param feature_id_column
#' @param ... additional arguments to plot_peptide_level function
#'
#' @return
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
plot_spike_ins <- function(data_df_long, sample_annotation,
                           spike_ins = 'BOVIN',
                           order_column = 'order',
                           measurement.col = 'Intensity',
                           batch_column = 'MS_batch',
                           sample_id_column = 'FullRunName',
                           feature_id_column = 'peptide_group_label', ...){
  spike_in_peptides = data_df_long %>%
    filter(grepl(spike_ins, ProteinName)) %>%
    pull(feature_id_column) %>% unique()
  gg = plot_peptide_level(spike_in_peptides, data_df_long, sample_annotation,
                          batch_column,feature_id_column, measurement.col, order_column,
                          title = 'Spike-in BOVINE protein peptides', ...)
  return(gg)
}

#' Plot Intensity for a few representative peptides for each step of the analysis including the fitting curve
#'
#' @param pep_name
#' @param data_df_all_steps
#' @param fit_df
#' @param sample_annotation
#' @param fit_value_var
#' @param batch_column
#' @param feature_id_column
#' @param measurement.col
#' @param sample_id_column
#' @param order_column
#' @param geom
#' @param color_by_batch
#' @param facet_by_batch
#' @param title
#' @param requant
#' @param theme
#' @param color_var
#'
#' @return
#' @export
#' @import ggplot2
#'
#' @examples
plot_with_fitting_curve <- function(pep_name, data_df_all_steps, fit_df, sample_annotation,
                                    fit_value_var = 'fit', fit_step = '3_loess_fit',
                                    batch_column = 'MS_batch',
                                    feature_id_column = 'peptide_group_label',
                                    measurement.col = 'Intensity',
                                    sample_id_column = 'FullRunName',
                                    order_column = NULL, geom = c('point', 'line'),
                                    color_by_batch = F, facet_by_batch = F,
                                    title = NULL, requant = NULL, theme = 'classic',
                                    color_var = 'fit'){
  if(length(pep_name) > 10){
    warning("Visualisation of individual features can be suboptimal,
            consider exploring no more than 5 features at a time")
  }
  gg = plot_peptide_level(pep_name, data_df_long = data_df_all_steps,
                          sample_annotation = sample_annotation,
                          sample_id_column = sample_id_column,
                          batch_column = batch_column,
                          feature_id_column = feature_id_column,
                          measurement.col = measurement.col,
                          order_column = order_column, title = title,
                          facet_by_batch = facet_by_batch)
  if(!("Step" %in% names(fit_df))){
    fit_df$Step = fit_step
    if (!(fit_step %in% data_df_all_steps$Step)){
      stop('specify step for which the curve fit should be shown')
    }
  }
  fit_df = fit_df %>%
    filter(rlang::UQ(as.name(feature_id_column)) %in% pep_name) %>%
    merge(sample_annotation, by = c(sample_id_column, batch_column))
  if(identical(color_by_batch, FALSE)){
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_var, x = order_column, group = batch_column),
                        color = 'red')+
      facet_grid(as.formula(paste(c(feature_id_column, "~", 'Step'), collapse = ' ')),
                 scales = 'free_y')
  } else {
    gg = gg + geom_line(data = fit_df,
                        aes_string(y = fit_value_var, x = order_column,
                                   group = batch_column, color = batch_column), size = 1.25)+
      facet_grid(as.formula(paste(c(feature_id_column, "~", 'Step'), collapse = ' ')),
                 scales = 'free_y')
    if(length(color_by_batch) == length(unique(fit_df[[batch_column]]))){
      gg = gg + scale_color_manual(values = color_by_batch)
    }
  }

  return(gg)
}
