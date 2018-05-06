#' Plot a single feature or several features each in its own facet
#'
#' single feature can also be spike-ins or iRT peptides (in DIA proteomics)
#'
#' @param pep_name name of the peptide for diagnostic profiling
#' @param df_long data frame where each row is a single feature in a single sample,
#' thus it has minimally, `sample_id_col`, `feature_id_column` and `measure_column`,
#' but usually also `m_score` (in OpenSWATH output result file)
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names)
#' 2) biological and 3) technical covariates (batches etc)
#' @param order_column sample injection order
#' @param sample_id_col name of the column in sample_annotation file,
#' where the filenames (colnames of the data matrix are found)
#' @param batch_column column in `sample_annotation` that should be used for batch comparison
#' @param measure_column if `df_long` is among the parameters, it is the column with
#' expression/abundance/intensity, otherwise, it is used internally for consistency
#' @param feature_id_column name of the column with feature/gene/peptide/protein ID
#' used with long format matrix (`df_long`).
#' Identical to wide format (`data_matrix`) to row name
#' @name feature_level_diagnostics


#' @name feature_level_diagnostics
#' @param geom whether to show the feature as points and/or connect by lines
#' @param color_by_batch (logical) whether to color points by batch
#' @param facet_by_batch (logical) whether to plot each batch in its own facet
#' @param title the string indicating the source of the peptides
#' @param requant if data frame: requant values;
#' if logical: whether to indicate requant values (requires 'requant' or 'm_score' column in `df_long`)
#' @param theme plot theme (default is 'classical'; other options not implemented)
#'
#' @return ggplot2 type plot of `measure_column` vs `order`,
#' faceted by peptide name and (optionally) by batch
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr %>%
#'
#' @examples
plot_peptide_level <- function(pep_name, df_long, sample_annotation,
                               order_column = NULL,
                               sample_id_col = 'FullRunName',
                               batch_column = 'MS_batch.final',
                               measure_column = 'Intensity',
                               feature_id_column = 'peptide_group_label',
                               geom = c('point', 'line'),
                               color_by_batch = F, facet_by_batch = F,
                               title = NULL, requant = NULL, theme = 'classic'){
  #TODO: suggest faceting by instrument
  #TODO: plot fit after LOESS



  plot_df = df_long %>%
    filter(rlang::UQ(as.name(feature_id_column)) %in% pep_name)
  if (!all(names(sample_annotation) %in% names(df_long))){
    sample_annotation = sample_annotation %>% arrange_(order_column)
    #remove annotation columns, except the sample_id_column
    common_cols = intersect(names(sample_annotation), names(plot_df))
    cols_to_remove = setdiff(common_cols, sample_id_col)
    plot_df = plot_df %>%
      select(-one_of(cols_to_remove))
    plot_df = plot_df %>%
      merge(sample_annotation, by = sample_id_col)
  }


  if(is.null(order_column)){
    warning("order column wasn't specified, putting row number as an order within a batch")
    plot_df = plot_df %>%
      group_by_at(vars(one_of(batch_column))) %>%
      mutate(order = row_number())
    order_column = 'order'
  }



  gg = ggplot(plot_df,
              aes_string(x = order_column, y = measure_column))
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
      filter(rlang::UQ(as.name(feature_id_column)) %in% pep_name)
    gg = gg + geom_point(data = data_requant,
                         aes_string(x = order_column, y = measure_column),
                         color = 'red', size = .3, shape = 8)
  }
  if (theme == 'classic'){
    gg = gg + theme_classic()
  }
  return(gg)
}

#' @name feature_level_diagnostics
#' @param spike_ins substring used to identify spike-in proteins in the column 'ProteinName'
#' @param ... additional arguments to plot_peptide_level function
#'
#' @export
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr %>%
#'
#' @examples
plot_spike_ins <- function(df_long, sample_annotation,
                           order_column = 'order',
                           spike_ins = 'BOVIN',
                           sample_id_column = 'FullRunName',
                           batch_column = 'MS_batch',
                           measure_column = 'Intensity',
                           feature_id_column = 'peptide_group_label',
                           title = 'Spike-in BOVINE protein peptides', ...){
  spike_in_peptides = df_long %>%
    filter(grepl(spike_ins, ProteinName)) %>%
    pull(feature_id_column) %>% unique()
  gg = plot_peptide_level(spike_in_peptides, df_long = df_long,
                          sample_annotation = sample_annotation,
                          order_column = order_column,
                          sample_id_col = sample_id_column,
                          batch_column = batch_column, measure_column = measure_column,
                          feature_id_column = feature_id_column,
                          title = title, ...)
  return(gg)
}

#' @name feature_level_diagnostics
#'
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#'
#' @examples
plot_iRTs <- function(df_long, sample_annotation,
                      order_column = NULL,
                      batch_column = 'MS_batch.final',
                      sample_id_col = 'FullRunName',
                      feature_id_column = 'peptide_group_label',
                      measure_column = 'Intensity',
                      title = 'iRT peptide profile', ...){
  iRT_peptides = df_long %>%
    filter(grepl('iRT', ProteinName)) %>%
    pull(feature_id_column)  %>% unique()
  gg = plot_peptide_level(iRT_peptides, df_long, sample_annotation,
                          order_column = order_column,
                          sample_id_col = sample_id_col,
                          batch_column = batch_column,
                          feature_id_column = feature_id_column,
                          measure_column = measure_column,
                          title = 'iRT_peptides', ...)
  return(gg)
}



#' Plot Intensity for a few representative peptides for each step of the analysis including the fitting curve
#'
#' @param pep_name name of the peptide for diagnostic profiling
#' @param data_df_all_steps data frame, similar to `df_long`,  where each row is a single feature in a single sample,
#' at a certain step of the analysis (minimally raw and after linear normalization)
#' thus it has minimally, `sample_id_col`, `feature_id_column`, `measure_column`  and `fit_step`,
#' but usually also `m_score` (in OpenSWATH output result file)
#' @param fit_df
#' @param fit_value_var
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names)
#' 2) biological and 3) technical covariates (batches etc)
#' @param order_column sample injection order
#' @param sample_id_col name of the column in sample_annotation file,
#' where the filenames (colnames of the data matrix are found)
#' @param batch_column column in `sample_annotation` that should be used for batch comparison
#' @param measure_column the column with expression/abundance/intensity.
#' @param feature_id_column name of the column with feature/gene/peptide/protein ID
#' used with long format matrix (`df_long`).
#' Identical to wide format (`data_matrix`) to row name
#' @param geom for the intensity (`measure_col`) profile:
#' @param color_by_batch
#' @param facet_by_batch
#' @param title
#' @param requant
#' @param theme
#'
#' @return `ggplot`-class plot with minimally two facets (before and after non-linear fit)
#' with `measure_column` (Intensity) vs `order_column` (injection order)
#' for selected peptides (specified in `pep_name`)
#' @export
#' @import ggplot2
#'
#' @examples
plot_with_fitting_curve <- function(pep_name, data_df_all_steps,
                                    sample_annotation,
                                    fit_df,
                                    fit_value_var = 'fit', fit_step = '3_loess_fit',
                                    order_column = NULL,
                                    sample_id_col = 'FullRunName',
                                    batch_column = 'MS_batch',
                                    measure_column = 'Intensity',
                                    feature_id_column = 'peptide_group_label',
                                    geom = c('point', 'line'),
                                    color_by_batch = F, facet_by_batch = F,
                                    title = NULL, requant = NULL,
                                    theme = 'classic'){
  #TODO: in non-linear fit, change "Intensity_normalized" to "Intensity", but rename "Intensity" as "Intensity_before_fit"
  if(length(pep_name) > 10){
    warning("Visualisation of individual features can be suboptimal,
            consider exploring no more than 5 features at a time")
  }
  gg = plot_peptide_level(pep_name, df_long = data_df_all_steps,
                          sample_annotation = sample_annotation,
                          order_column = order_column,
                          sample_id_col = sample_id_col,
                          batch_column = batch_column,
                          measure_column = measure_column,
                          feature_id_column = feature_id_column,
                          title = title,
                          facet_by_batch = facet_by_batch)
  if(!("Step" %in% names(fit_df))){
    fit_df$Step = fit_step
    if (!(fit_step %in% data_df_all_steps$Step)){
      stop('specify step for which the curve fit should be shown')
    }
  }
  fit_df = fit_df %>%
    filter(rlang::UQ(as.name(feature_id_column)) %in% pep_name) %>%
    merge(sample_annotation, by = c(sample_id_col, batch_column))
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
