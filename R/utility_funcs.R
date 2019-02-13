check_sample_consitency <- function(sample_annotation, sample_id_col, df_long) {
  if (!is.null(sample_annotation)){
    if (!(sample_id_col %in% names(sample_annotation))){
      warning('Sample ID column is not defined in sample annotation, sample annotation 
              cannot be used for initial assessment of batch effects')
    }
    if(!setequal(unique(sample_annotation[[sample_id_col]]), unique(df_long[[sample_id_col]]))){
      warning('Sample IDs in sample annotation not consistent with samples in input data, 
              will merge, using intersecting Sample IDs only')
    }
    df_long = df_long %>% merge(sample_annotation,
                                by = sample_id_col)
  } else {
    warning('Sample annotation is not provided, only the basic sample boxplots will be plotted')
  }
}

define_sample_order <- function(order_col, sample_annotation, facet_col, batch_col, 
                                df_long, sample_id_col, color_by_batch) {
  if (!is.null(order_col)){
    if (!is.null(sample_annotation) & !(order_col %in% names(sample_annotation))){
      if(!is.null(facet_col)){
        warning("order column not found in sample annotation, assuming that order in 
                sample annotation corresponds to sample running order, specific for each instrument")
      } else {
        warning('order column not found in sample annotation, 
                assuming the order of sample IDs corresponds to running order')
      }
    } else {
      warning('sample annotation is not defined, 
                taking order of files in the data matrix instead')
    }
  } else {
    if (!is.null(batch_col)){
      warning("Order column is NULL, assuming order is not introducing unwanted 
            association between the samples, plotting samples in order of batch factor")
      df_long[[batch_col]] = as.factor(df_long[[batch_col]])
    } else {
      warning("Order column is NULL, assuming order is not introducing unwanted 
            association between the samples")
    }
  }
  
  if (order_col == sample_id_col & color_by_batch & (batch_col %in% names(sample_annotation))){
    if (!is.null(sample_annotation)){
      warning('order column is identical to sample ID column, 
                assuming order of samples in the annnotation corresponds to the sample running order')
      order_col = 'sample_order'
      df_long[[order_col]] = match(df_long[[sample_id_col]],
                                  sample_annotation[[sample_id_col]])
    } else {
      warning('order column is identical to sample ID column, and sample annotation is not defined,
                assuming order of samples in the intensity table corresponds to the sample running order')
      order_col = 'sample_order'
      df_long[[order_col]] = match(df_long[[sample_id_col]],
                                  unique(df_long[[sample_id_col]]))
    }
  }
  
  if (is.null(order_col) | !(order_col %in% names(df_long))){
    order_col = sample_id_col
    df_long[[order_col]] = factor(df_long[[order_col]], levels = df_long[[order_col]])
  }
  
  #infer the order within facets
  if(!is.null(facet_col) & is.numeric(df_long[[order_col]])){
    df_long = df_long %>% 
      group_by_at(vars(one_of(facet_col))) %>% 
      mutate(order_per_instrument = dense_rank(UQ(sym(order_col))))
    order_col = 'order_per_instrument'
  }
  
  return(list(order_col = order_col,
              df_ave = df_long))
}

add_vertical_batch_borders <- function(order_col, batch_col, vline_color, facet_col, 
                                       df_long, gg) {
  if(!is.null(order_col) & !is.null(batch_col) & !is.null(vline_color)){
    #define the batch tipping points (positions of vertical lines)
    if (!is.null(facet_col)){
      order_vars <- c(facet_col, order_col)
      batch_vars = c(facet_col, batch_col)
      batch.tipping.points = df_long %>%
        arrange(!!!syms(order_vars))%>%
        group_by(!!!syms(batch_vars)) %>%
        summarise(batch_size = n()) %>%
        group_by(!!sym(facet_col)) %>%
        mutate(tipping.points = cumsum(batch_size))%>%
        mutate(tipping.points = tipping.points+.5)
    } else {
      batch.tipping.points = df_long %>%
        arrange(!!sym(order_col))%>%
        group_by(!!sym(batch_col)) %>%
        summarise(batch_size = n()) %>%
        mutate(tipping.points = cumsum(batch_size))%>%
        mutate(tipping.points = tipping.points+.5)
    }
    gg = gg + geom_vline(data = batch.tipping.points,
                         aes(xintercept = tipping.points),
                         color = vline_color, linetype = 'dashed')
  }
  return(gg)
}