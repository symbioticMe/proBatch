check_sample_consistency <- function(sample_annotation, sample_id_col, df_long) {
  if (!is.null(sample_annotation)){
    if (!(sample_id_col %in% names(sample_annotation))){
      warning('Sample ID column is not defined in sample annotation, sample annotation 
              cannot be used for initial assessment of batch effects')
    }
    if(!setequal(unique(sample_annotation[[sample_id_col]]), unique(df_long[[sample_id_col]]))){
      warning('Sample IDs in sample annotation not consistent with samples in input data, 
              will merge, using intersecting Sample IDs only')
      #TODO: expand the warnings for more specific cases: 1) sample annotation has samples not represented in data matrix; 2) dm has samples not in annotation;
      #TODO: Break the merge if 1) sample annotation has duplicated samples; 2) dm has duplicated samples
    }
    df_long = df_long %>% 
      inner_join(sample_annotation, by = sample_id_col) %>%
      as.data.frame()
  } else {
    warning('Sample annotation is not provided, only the basic plot will be visualized')
  }
  return(df_long)
}

#' Defining sample order internally 
#'
#' @inheritParams proBatch
#'
#' @return   list of two items: \code{order_col} new name and new \code{df_long}
#' @export
#'
#' @examples 
#' sample_order = define_sample_order(order_col, sample_annotation, 
#' facet_col, batch_col, df_long, sample_id_col, color_by_batch)
#' new_order_col = sample_order$order_col
#' df_long = sample_order$df_long
#' 
#' @seealso \link{plot_sample_mean_or_boxplot}, \link{feature_level_diagnostics}
#' 
#' @name define_sample_order
define_sample_order <- function(order_col, sample_annotation, facet_col, batch_col, 
                                df_long, sample_id_col, color_by_batch) {
  if (!is.null(order_col)){
    if (!is.null(sample_annotation)){
      if (!(order_col %in% names(sample_annotation))){
        if(!is.null(facet_col)){
          warning(sprintf('column  %s is not found in sample annotation, assuming that order in 
                sample annotation corresponds to sample running order, specific for each instrument', order_col))
        } else {
          warning(sprintf('column %s is not defined in sample annotation, 
                assuming the order of sample IDs corresponds to running order', order_col))
        }
      } else {
        if (order_col != sample_id_col){
          return(list(order_col = order_col,
                      df_long = df_long))
        }
      }
    } else {
      if (!(order_col %in% names(df_long))){
        if(!is.null(facet_col)){
          warning(sprintf('column  %s for order is not found in data frame, 
                          taking order of files in the data matrix instead, specific for each instrument', order_col))
        } else {
          warning(sprintf('column %s is not defined in data frame, 
                taking order of files in the data matrix instead', order_col))
        }
      } else {
        return(list(order_col = order_col,
                    df_long = df_long))
      }
    } 
  } else {
    if (!is.null(batch_col) && (batch_col %in% names(df_long))){
      warning("Order column is NULL, assuming order is not introducing unwanted 
            association between the samples, plotting samples in order of batch factor")
      df_long[[batch_col]] = as.factor(df_long[[batch_col]])
    } else {
      warning("Order column is NULL, assuming order is not introducing unwanted 
            association between the samples")
    }
  }
  
  if (!is.null(order_col) && order_col == sample_id_col){
    if (!is.null(sample_annotation)){
      order_col = 'sample_order'
      if(color_by_batch && (!is.null(batch_col) && (batch_col %in% names(sample_annotation)))){
        warning("order column is identical to sample ID column and coloring by batch is required,
                ordering the samples by batch rather than by sample order in annotation")
        df_long = df_long %>% arrange(!!sym(c(batch_col))) %>%
          mutate(UQ(sym(order_col)) :=  factor(UQ(sym(sample_id_col)),
                                                      levels = unique(UQ(sym(sample_id_col)))))
        #df_long[[order_col]] = reorder(as.character(df_long[[sample_id_col]]), df_long[[batch_col]])
      } else {
        warning('order column is identical to sample ID column, 
                assuming order of samples in the annnotation corresponds to the sample running order')
        df_long[[order_col]] = match(df_long[[sample_id_col]],
                                     sample_annotation[[sample_id_col]])
      }
    } else {
      warning('order column is identical to sample ID column, and sample annotation is not defined,
                assuming order of samples in the intensity table corresponds to the sample running order')
      order_col = 'sample_order'
      df_long[[order_col]] = match(df_long[[sample_id_col]],
                                   unique(df_long[[sample_id_col]]))
    }
  }
  
  if (is.null(order_col) || !(order_col %in% names(df_long))){
    order_col = sample_id_col
    if(color_by_batch & (!is.null(batch_col) && (batch_col %in% names(sample_annotation)))){
      warning("order column is not defined and coloring by batch is required,
                ordering the samples by batch")
      df_long = df_long %>% arrange(!!sym(c(batch_col))) %>%
        mutate(UQ(sym(order_col)) :=  factor(UQ(sym(order_col)),
                                                    levels = unique(UQ(sym(order_col)))))
      #df_long[[order_col]] = reorder(as.character(df_long[[order_col]]), df_long[[batch_col]])
    } else {
      df_long[[order_col]] = factor(df_long[[order_col]], levels = unique(df_long[[order_col]]))
    }
  }
  
  #infer the order within facets
  if(!is.null(facet_col) && is.numeric(df_long[[order_col]])){
    df_long = df_long %>% 
      group_by_at(vars(one_of(facet_col))) %>% 
      mutate(order_per_instrument = dense_rank(UQ(sym(order_col))))
    order_col = 'order_per_instrument'
  }
  
  return(list(order_col = order_col,
              df_long = df_long))
}

add_vertical_batch_borders <- function(order_col, sample_id_col, batch_col, vline_color, facet_col, 
                                       sample_annotation, gg) {
  if(!is.null(order_col) & (order_col != sample_id_col) & 
     !(is.character(sample_annotation[[order_col]]) || is.factor(sample_annotation[[order_col]]))&
     !is.null(batch_col) & !is.null(vline_color)){
    #define the batch tipping points (positions of vertical lines)
    if (!is.null(facet_col)){
      sample_annotation = sample_annotation %>%
        select(one_of(c(order_col, sample_id_col, batch_col, facet_col))) %>%
        distinct()
      order_vars <- c(facet_col, order_col)
      batch_vars = c(facet_col, batch_col)
      min_order_values = sample_annotation %>%
        group_by(!!sym(facet_col)) %>%
        summarise(min_order_value = min(!!sym(order_col))-1)
      batch.tipping.points = sample_annotation %>%
        arrange(!!!syms(order_vars))%>%
        group_by(!!!syms(batch_vars)) %>%
        summarise(batch_size = n()) %>%
        ungroup() %>%
        merge(min_order_values, by = facet_col) %>%
        group_by(!!sym(facet_col)) %>%
        mutate(tipping.points = cumsum(batch_size))%>%
        mutate(tipping.points = tipping.points+.5 + min_order_value)
    } else {
      sample_annotation = sample_annotation %>%
        select(one_of(c(order_col, sample_id_col, batch_col))) %>%
        distinct()
      min_order_val = min(sample_annotation[[order_col]]) - 1
      batch.tipping.points = sample_annotation %>%
        arrange(!!sym(order_col))%>%
        group_by(!!sym(batch_col)) %>%
        summarise(batch_size = n()) %>%
        ungroup() %>%
        mutate() %>% #enables adequate plotting for batches, where order doesn't start from one
        mutate(tipping.points = cumsum(batch_size))%>%
        mutate(tipping.points = tipping.points+.5 + min_order_val)
    }
    gg = gg + geom_vline(data = batch.tipping.points,
                         aes(xintercept = tipping.points),
                         color = vline_color, linetype = 'dashed')
  }
  return(gg)
}