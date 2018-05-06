map_factors_to_colors <- function(annotation_df_factors){

  #calculate number of colors to create
  nlev_covariate = mapply(nlevels, annotation_df_factors)
  n_colors_total = sum(nlev_covariate)

  if (any(nlev_covariate > 20)){
    warning('Some colors will be hard to distinguish\n')
  }
  if(any(nlev_covariate > 50)){
    warning('Too many colors, consider merging some covariate values for better visualisation\n')
  }

  colors = WGCNA::standardColors(n_colors_total)
  start_indxs = c(1, 1 + cumsum(nlev_covariate[-length(nlev_covariate)]))
  end_indx = cumsum(nlev_covariate)
  ann_colors_covariate = lapply(1:length(nlev_covariate),
                                function(i) colors[start_indxs[i]:end_indx[i]])

  names(ann_colors_covariate) = names(nlev_covariate)
  ann_colors_covariate <- Map(setNames, ann_colors_covariate, lapply(annotation_df_factors, levels))
  return(ann_colors_covariate)
}

map_numbers_to_colors <- function(annotation_df_numbers, palette_type = 'brewer',
                                  granularity = 10){

  n_colors_to_create <- ncol(annotation_df_numbers)
  if ((n_colors_to_create > 4 & palette_type == 'viridis')){
    warning('Too many colors for viridis palette, switching to Brewer palettes')
  }
  if ((n_colors_to_create > 18)){
    stop('Not enough color paletters to visualize the annotation')
  }

  if (length(granularity) > 1 & length(granularity) != ncol(annotation_df_numbers)){
    warning('Either specify universal color granularity or granularity for each of numeric vectors!')
    granularity = granularity[1]
  }

  if (length(granularity) == 1) {
    ann_col_covariate = lapply(1:ncol(annotation_df_numbers),
                               function(i) generate_colors_for_numeric(annotation_df_numbers[, i],
                                                                       i = i,
                                                                       palette_type = palette_type,
                                                                       granularity = granularity))
  } else {
    ann_col_covariate = lapply(1:ncol(annotation_df_numbers),
                               function(i) generate_colors_for_numeric(annotation_df_numbers[, i],
                                                                       i = i,
                                                                       palette_type = palette_type,
                                                                       granularity = granularity[i]))
  }

  color_list = lapply(ann_col_covariate, function(item) item$color_vector)
  names(color_list) = names(annotation_df_numbers)

  new_sample_annotation = data.frame(lapply(ann_col_covariate, function(item) item$new_annotation))
  names(new_sample_annotation) = names(annotation_df_numbers)

  return(list(color_list = color_list,
              new_sample_annotation = new_sample_annotation))

}

#' generate a list of colors for the dataframe with all columns numeric (or date)
#'
#' @param num_col
#' @param palette_type
#' @param column_to_log
#' @param i
#'
#' @return
#' @importFrom scales brewer_pal
#' @importFrom lubridate is.POSIXct
#'
#' @examples
generate_colors_for_numeric <- function(num_col, palette_type = 'brewer',
                                        column_to_log = F, i = 1, granularity = 10){
  color_for_column = switch(palette_type,
                            brewer = brewer_pal(type = "div", i)(5)[1:5],
                            viridis = viridis::viridis_pal(option = LETTERS[5-i])(5))

  non_numeric_values = NULL
  if(is.factor(num_col)){
    num_col_temp = as.character(num_col)
    if(all(is.na(as.numeric(num_col_temp)))){
      warning('This column is not a numeric')
    }
    num_col = as.numeric(num_col_temp)
    non_numeric_values = unique(num_col_temp[is.na(num_col)])
    n_non_numeric = length(non_numeric_values)
    factor_colors = scales::brewer_pal(palette = 'Set1')(n_non_numeric)
    names(factor_colors) = non_numeric_values
  }

  if (is.numeric(num_col)){
    num_vec = cut(num_col, breaks = granularity)
  } else if (is.POSIXct(num_col)) {
    interval = (max(num_col, na.rm = T) - min(num_col, na.rm = T))/granularity
    if(any(is.na(num_col))){
      warning('NAs in the numeric vector')
      }
    interval_char = paste(as.integer(interval), units(interval), sep = ' ')
    num_vec = cut(num_col, breaks = interval_char)
  }

  #color_to_plot = colorRampPalette(color_for_column)(nlevels(num_vec))[1:nlevels(num_vec)]
  color_to_plot = colorRampPalette(color_for_column)(nlevels(num_vec))[1:nlevels(num_vec)]
  names(color_to_plot) = levels(num_vec)

  if(!is.null(non_numeric_values)){
    #merge the factors back together
    num_vec = as.character(num_vec)
    num_vec[is.na(num_vec)] = num_col_temp[is.na(num_vec)]
    num_vec = as.factor(num_vec)
    color_to_plot = c(color_to_plot, factor_colors)
  }

  return(list(color_vector = color_to_plot,
              new_annotation = num_vec))
}

check_rare_levels <- function(col){
  tb_col = table(col)
  freq_table = table(tb_col)/length(tb_col)
  is_rare = as.character(1) %in% names(freq_table) & freq_table[as.character(1)] > .5
  return(is_rare)
}

merge_rare_levels <- function(col) {
  is_factor_col = is.factor(col)
  tb_col = table(col)
  if(is_factor_col) col = as.character(col)
  col[col %in% names(tb_col)[tb_col == 1]] = 'other'
  if(is_factor_col) col = as.factor(col)
  return(col)
}

#' convert the sample annotation data frame to list of colors
#' the list is named as columns included to use in potting functions
#'
#' @param columns_for_plotting
#' @param sample_id_column
#' @param factor_columns
#' @param not_factor_columns
#' @param rare_categories_to_other
#' @param numerics_to_log
#' @param granularity number of colors to map to the number vector (equally spaced between minimum and maximum)
#' @param sample_annotation
#'
#' @return list of colors
#' @export
#' @import dplyr
#' @import WGCNA
#' @importFrom magrittr %>%
#'
#' @examples
sample_annotation_to_colors <- function(sample_annotation,
                                        columns_for_plotting = NULL,
                                        sample_id_column = NULL,
                                        factor_columns = NULL,
                                        not_factor_columns = NULL,
                                        rare_categories_to_other = T,
                                        numerics_to_log = F,
                                        numeric_palette_type = 'brewer',
                                        granularity = 10){
  rownames_ann = as.character(sample_annotation[[sample_id_column]])
  if(!is.null(columns_for_plotting)){
    sample_annotation = sample_annotation %>%
      select(one_of(setdiff(columns_for_plotting, sample_id_column)))
  }

  factor_or_not <- intersect(factor_columns, not_factor_columns)
  if(length(factor_or_not) > 1){
    stop(sprintf('Columns %s are defined as both factors and not factors',
                 paste(factor_or_not, collapse = ', ')))
  }

  factor_like_columns = names(sample_annotation)[sapply(sample_annotation,
                               function(col)
                                 is.factor(col) | is.character(col))]
  if(!is.null(factor_columns)){
    factor_columns = union(factor_columns, factor_like_columns)
    }
  if (!is.null(not_factor_columns)){
    factor_columns = setdiff(factor_columns, not_factor_columns)
  }

  #TODO: check if this is absolutely required (convertion to factors)
  list_of_col_for_factors = list()
  if(!is.null(factor_columns)){
    factor_df = sample_annotation %>%
      select(one_of(factor_columns)) %>%
      mutate_at(names(.), funs(as.factor))

    if (rare_categories_to_other){
      factor_df = factor_df %>%
        mutate_if(check_rare_levels, funs(merge_rare_levels(.)))
    }

    #generate factor mappings
    list_of_col_for_factors = map_factors_to_colors(factor_df)
  }


  #generate color mappings of numeric variables
  non_factor_cols = setdiff(names(sample_annotation), factor_columns)
  #TODO: if numerics_to_log is a character vector of column names, convert corresponding annotation colors to log scale
  list_of_col_for_numeric = list()
  if(!is.null(non_factor_cols)){
    numeric_df = sample_annotation %>%
      select(one_of(non_factor_cols))
    map_of_colors_to_num_vec = map_numbers_to_colors(numeric_df,
                                                     palette_type = numeric_palette_type,
                                                     granularity = granularity)
    list_of_col_for_numeric = map_of_colors_to_num_vec$color_list
    numeric_factor_df = map_of_colors_to_num_vec$new_sample_annotation
    sample_annotation = cbind(factor_df, numeric_factor_df)
  }

  #join two lists of colors
  list_of_colors = c(list_of_col_for_factors, list_of_col_for_numeric)
  rownames(sample_annotation) = rownames_ann
  return(list(list_of_colors = list_of_colors,
         #TODO: df_of_colors = df_of_colors,
         sample_annotation = sample_annotation))
}

#' Turn color list to df (some plotting functions require the latter)
#'
#' @param color_list
#' @param sample_annotation
#'
#' @return
#' @export
#'
#' @examples
color_list_to_df <- function(color_list, sample_annotation){
  list_df = lapply(names(sample_annotation), function(col_name){
    col_values = sample_annotation[[col_name]]
    col_colors = color_list[[col_name]][col_values]
  })
  names(list_df) = names(sample_annotation)
  color_df = as.data.frame(do.call(cbind, list_df))
  rownames(color_df) = rownames(sample_annotation)
  return(color_df)
}