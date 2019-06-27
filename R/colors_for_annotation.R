map_factors_to_colors <- function(annotation_df_factors) {
  #calculate number of colors to create
  col_class = vapply(annotation_df_factors, 
                     FUN = function(x) paste(class(x), collapse = "; "),
                     FUN.VALUE = character(1))
  non_factor_cols = names(annotation_df_factors)[col_class != 'factor']
  wrong_classes = col_class[col_class != 'factor']
  if (length(non_factor_cols) > 0){
    col_string = paste(non_factor_cols, collapse = ', ')
    col_classes = paste(wrong_classes, collapse = ', ')
    warning(sprintf('Columns %s are non factors, but %s, they will be converted 
                    to factors for mapping to colors', col_string, col_classes))
    annotation_df_factors = annotation_df_factors %>% 
      mutate_if(negate(is.factor), as.factor)
  }
  
  nlev_covariate = mapply(nlevels, annotation_df_factors)
  n_colors_total = sum(nlev_covariate)
  
  if (any(nlev_covariate > 20)) {
    warning('Some colors will be hard to distinguish\n')
  }
  if (any(nlev_covariate > 50)) {
    warning('Too many colors, consider merging some covariate 
                values for better visualisation\n')
  }

  colors = standardColors(n_colors_total)
  start_indxs = c(1, 1 + cumsum(nlev_covariate[-length(nlev_covariate)]))
  end_indx = cumsum(nlev_covariate)
  ann_colors_covariate = lapply(seq_len(length(nlev_covariate)),
                                function(i)
                                  colors[start_indxs[i]:end_indx[i]])

  names(ann_colors_covariate) = names(nlev_covariate)
  ann_colors_covariate <-
      Map(setNames,
          ann_colors_covariate,
          lapply(annotation_df_factors, levels))
  change_other_to_grey <- function(x) { 
    if ('other' %in% names(x)){
      x['other'] <- 'grey'
    }
      return(x)
  }
  ann_colors_covariate = lapply(ann_colors_covariate, change_other_to_grey)
  return(ann_colors_covariate)
}

map_numbers_to_colors <-
  function(annotation_df_numbers,
           palette_type = 'brewer',
           granularity = 10) {
    n_colors_to_create <- ncol(annotation_df_numbers)
    if ((n_colors_to_create > 4 & palette_type == 'viridis')) {
      warning('Too many colors for viridis palette, 
                    switching to Brewer palettes')
    }
    if ((n_colors_to_create > 18)) {
      stop('Not enough color paletters to visualize the annotation')
    }
    
    if (length(granularity) > 1 &
        length(granularity) != ncol(annotation_df_numbers)) {
      warning(
        'Either specify universal color granularity or 
                granularity for each of numeric vectors!'
      )
      granularity = granularity[1]
    }
    
    if (length(granularity) == 1) {
      granularity <- rep(granularity, ncol(annotation_df_numbers))}
    
    ann_col_covariate <- lapply(seq_len(ncol(annotation_df_numbers)),
                                function(i)
                                  generate_colors_for_numeric(
                                    annotation_df_numbers[, i],
                                    i = i,
                                    palette_type = palette_type,
                                    granularity = granularity[i]
                                  ))
    
    color_list = lapply(ann_col_covariate, function(item)
      item$color_vector)
    names(color_list) = names(annotation_df_numbers)
    new_sample_annotation = data.frame(lapply(
      ann_col_covariate, function(item)item$new_annotation))
    names(new_sample_annotation) = names(annotation_df_numbers)
    
    return(list(color_list = color_list,
                new_sample_annotation = new_sample_annotation))
    
  }


#' Generates color list
#'
#' Generates a list of colors for a vector of numeric, POSIXct (i.e. the
#' (signed) number of seconds since the beginning of 1970 , or factors
#'
#' @param num_col a vector of type numeric of factor to generate colors for
#' @param palette_type 'brewer' or 'viridis'
#' @param i if \code{palette_type} is 'brewer' the palette argument to
#'   \code{brewer_pal}. If \code{palette_type} is 'viridis' the option argument
#'   to \code{viridis_pal}
#' @param granularity breaks to use when generating colors for \code{num_col}
#'
#' @return list, containing the following items:
#' \enumerate{
#' \item \code{color_vector} - string-like vector of colors
#' \item \code{new_annotation} - factor representation of numeric vector 
#' (factor with number of levels equal to "granularity")
#' }
#' @keywords internal
#' 
generate_colors_for_numeric <- function(num_col,
                                        palette_type = 'brewer',
                                        i = 1,
                                        granularity = 10) {
  if ((palette_type == 'viridis') && (i > 4 || i < 1)) {
    warning('When using viridis palette i 
                must be >= 1 and <= 4. Setting it to 1.')
    i = 1
  }
  
  color_for_column = switch(
    palette_type,
    brewer = brewer_pal(type = "div", i)(5)[seq_len(5)],
    viridis = viridis_pal(option = LETTERS[5 - i])(5)
  )
  
  non_numeric_values = NULL
  if (is.factor(num_col)) {
    num_col_temp = as.character(num_col)
    if (all(is.na(as.numeric(num_col_temp)))) {
      warning('This column is not a numeric')
    }
    num_col = as.numeric(num_col_temp)
    non_numeric_values = unique(num_col_temp[is.na(num_col)])
    n_non_numeric = length(non_numeric_values)
    factor_colors = brewer_pal(palette = 'Set1')(n_non_numeric)
    names(factor_colors) = non_numeric_values
  }
  
  if (is.numeric(num_col)) {
    num_vec = cut(num_col, breaks = granularity)
  } else if (is.POSIXct(num_col)) {
    interval = (max(num_col, na.rm = TRUE) - 
                  min(num_col, na.rm = TRUE)) / granularity
    if (any(is.na(num_col))) {
      warning('NAs in the numeric vector')
    }
    interval_char = paste(as.integer(interval), units(interval), sep = ' ')
    num_vec = cut(num_col, breaks = interval_char)
  }
  
  color_to_plot = colorRampPalette(color_for_column)(
    nlevels(num_vec))[seq_len(nlevels(num_vec))]
  names(color_to_plot) = levels(num_vec)
  
  if (!is.null(non_numeric_values)) {
    num_vec = as.character(num_vec)
    num_vec[is.na(num_vec)] = num_col_temp[is.na(num_vec)]
    num_vec = as.factor(num_vec)
    color_to_plot = c(color_to_plot, factor_colors)
  }
  
  return(list(color_vector = color_to_plot,
              new_annotation = num_vec))
}

check_rare_levels <- function(column) {
  tb_col = table(column)
  freq_table = table(tb_col) / length(tb_col)
  is_rare = as.character(1) %in% names(freq_table) &
    freq_table[as.character(1)] > .5
  return(is_rare)
}


#' Replaces rare levels with other
#'
#' Replaces levels with a maximal occurrence of 1 with \code{other}
#'
#' @keywords internal
#'
#' @return column with rare occurrences replaced by other
merge_rare_levels <- function(column) {
  is_factor_col = is.factor(column)
  tb_col = table(column)
  if (is_factor_col)
    column = as.character(column)
  column[column %in% names(tb_col)[tb_col == 1]] = 'other'
  if (is_factor_col)
    column = as.factor(column)
  return(column)
}


#' Generate colors for sample annotation
#'
#' Convert the sample annotation data frame to list of colors
#' the list is named as columns included to use in plotting functions
#'
#' @inheritParams proBatch
#' @param factor_columns columns of \code{sample_annotation} to be 
#' treated as factors. Sometimes categorical variables are depicted as integers 
#' (e.g. in column "Batch", values are 1, 2 and 3), specification here allows to
#' map them correctly to qualitative palettes.
#' @param date_columns columns of `POSIXct` class. These columns are mapped to 
#' continuous palettes, as numeric columns, but require special treatment.
#' @param numeric_columns columns of \code{sample_annotation} to be 
#' treated as continuous numeric values. 
#' @param rare_categories_to_other if \code{True} rare categories 
#' will be merged into the value \code{"other"}
#' @param granularity number of colors to map to the number 
#' vector (equally spaced between minimum and maximum)
#' @param numeric_palette_type palette to be used for 
#' numeric values coloring
#'
#' @return list of three items: \enumerate{
#'   \item list of colors; 
#'   \item data frame of colors; 
#'   \item new sample annotation (e.g. rare factor levels merged into "other")
#'   }
#' 
#' @examples 
#' color_scheme <- sample_annotation_to_colors (example_sample_annotation, 
#' factor_columns = c('MS_batch','EarTag', "Strain", 
#' "Diet", "digestion_batch", "Sex"),
#' date_columns = 'DateTime',
#' numeric_columns = c('order'))
#' @export
#' 
#' @name sample_annotation_to_colors
sample_annotation_to_colors <- function(sample_annotation,
                                        sample_id_col = 'FullRunName',
                                        factor_columns = c('MS_batch','EarTag', 
                                                           'digestion_batch',
                                                           "Strain", "Diet"),
                                        date_columns = 'DateTime',
                                        numeric_columns = 'order',
                                        rare_categories_to_other = TRUE,
                                        numeric_palette_type = 'brewer',
                                        granularity = 10) {
  
  sample_annotation = as.data.frame(sample_annotation)
  message('converting columns to corresponding classes 
          (factor, POSIXct date, numeric)')
  sample_annotation <- sample_annotation %>%
    mutate_at(vars(factor_columns), as.factor) %>%
    mutate_at(vars(numeric_columns), as.numeric) %>%
    mutate_at(vars(date_columns), as.POSIXct, date_columns)
  
  columns_for_color_mapping = union(factor_columns, 
                                    union(numeric_columns, date_columns))
  undefined_cols <- setdiff(names(sample_annotation), 
                            c(columns_for_color_mapping, sample_id_col))
  if (length(undefined_cols) > 0){
    warning(paste(c('The following columns will not be mapped to colors:',
                    undefined_cols, 'if these have to be mapped, please assign 
                    them to factor, date or numeric'), collapse = ' '))
  }
  
  rownames_ann = as.character(sample_annotation[[sample_id_col]])
  sample_annotation = sample_annotation %>%
    select(one_of(columns_for_color_mapping))
  
  if (!is.null(date_columns) || !is.null(numeric_columns)) {
    factor_columns = setdiff(factor_columns, c(date_columns, numeric_columns))
    column_intersection <- intersect(factor_columns, 
                                     union(date_columns, numeric_columns))
    if (length(column_intersection) > 0) {
      warning(paste(c('The following columns are repeatedly listed among factors
                      and numeric-like variables:', column_intersection, 
                      '; they will be excluded from factors and mapped to 
                      continuous palettes'), 
                    collapse = ' '))
    }
  }
  
  list_of_col_for_factors = list()
  if (!is.null(factor_columns)) {
    factor_df = sample_annotation %>%
      select(one_of(factor_columns))
    
    if (rare_categories_to_other) {
      factor_df = factor_df %>%
        mutate_if(check_rare_levels, funs(merge_rare_levels(.)))
    }
    
    list_of_col_for_factors = map_factors_to_colors(factor_df)
  }
  
  non_factor_cols = setdiff(names(sample_annotation), factor_columns)
  
  list_of_col_for_numeric = list()
  if (!is.null(non_factor_cols) &
      !identical(non_factor_cols, character(0))) {
    numeric_df = sample_annotation %>%
      select(one_of(non_factor_cols))
    map_of_colors_to_num_vec = map_numbers_to_colors(
      numeric_df,palette_type = numeric_palette_type,
      granularity = granularity)
    list_of_col_for_numeric = map_of_colors_to_num_vec$color_list
    numeric_factor_df = map_of_colors_to_num_vec$new_sample_annotation
    sample_annotation = cbind(factor_df, numeric_factor_df)
  } else {
      sample_annotation = factor_df
  }
  
  list_of_colors = c(list_of_col_for_factors, list_of_col_for_numeric)
  rownames(sample_annotation) = rownames_ann
  
  
  
  return(list(list_of_colors = list_of_colors,
              sample_annotation = sample_annotation))
}


#' Color list to data frame
#'
#' Turn color list to df (some plotting functions require the latter)
#'
#' @param color_list list of colors
#' @param sample_annotation factor-based configuration 
#' of the sample annotation
#' 
#' @return a data frame representation of the input color list
#' 
#' @keywords internal
#'
color_list_to_df <- function(color_list, sample_annotation) {
    list_df = lapply(names(sample_annotation), function(col_name) {
        col_values = sample_annotation[[col_name]]
        col_colors = color_list[[col_name]][col_values]
    })
    names(list_df) = names(sample_annotation)
    color_df = as.data.frame(do.call(cbind, list_df))
    rownames(color_df) = rownames(sample_annotation)
    return(color_df)
}

color_points_by_batch <- function(color_by_batch, batch_col, gg, color_scheme, 
                                  sample_annotation) {
  if(color_by_batch & !is.null(batch_col)){
    gg = gg + aes(color = !!sym(batch_col))
    
    n_batches <- length(unique(sample_annotation[[batch_col]]))
    if (n_batches == length(color_scheme)){
      
    }
    #Define the color scheme on the fly
    if(length(color_scheme) == 1 && color_scheme == 'brewer'){
      
      if (n_batches <= 9){
        gg = gg + scale_color_brewer(palette = 'Set1')
      } else {
        if (n_batches <= 12){
          gg = gg + scale_color_brewer(palette = 'Set3')
        } else {
          warning(sprintf('brewer palettes have maximally 12 colors, 
                           %s batches are specified,
                          consider defining color scheme with 
                          sample_annotation_to_colors function', 
                          n_batches))
        }
      }
    } else {
      #color vector provided by "sample_annotation_to_colors"
      gg = gg + scale_color_manual(values = color_scheme)
    }
  } else {
    if (is.null(batch_col)){
      stop('Coloring column not defined, please define the color column!')
    }
  }
  return(gg)
}

color_fill_boxes_by_batch <- function(color_by_batch, batch_col, gg, 
                                      color_scheme, df_long) {
  if(color_by_batch & !is.null(batch_col)){
    gg = gg + aes(fill = !!sym(batch_col))
    
    
    if(length(color_scheme) == 1 && color_scheme == 'brewer'){
      ##TODO: 1) check the type of the factor; 2) invoke discrete/continuous coloring
      n_batches <- length(unique(df_long[[batch_col]]))
      if(n_batches < 9){
        gg = gg + scale_fill_brewer(palette = 'Set1')
        
      } else {
        if (n_batches <= 12){
          gg = gg + scale_fill_brewer(palette = 'Set3')
        } else {
          warning(sprintf('brewer palettes have maximally 12 colors, 
                           %s batches are specified,
                          consider defining color scheme with 
                          sample_annotation_to_colors function', 
                          n_batches))
        }
        }
      } else{
        gg = gg + scale_fill_manual(values = color_scheme)
    }
  }
  return(gg)
}

color_discrete <- function(color_scheme, batch_col, n_batches, fill_or_color, gg) {
  
  if(fill_or_color == 'color'){
    gg = gg + aes(color = as.factor(!!sym(batch_col)))
  } else {
    if(fill_or_color == 'fill'){
      gg = gg + aes(fill = as.factor(!!sym(batch_col)))
    }
  }
  
  #Define the color scheme on the fly
  if(length(color_scheme) == 1 && color_scheme == 'brewer'){
    if (n_batches <= 9){
      if(fill_or_color == 'color'){
        gg = gg + scale_color_brewer(palette = 'Set1')
      } else {
        if(fill_or_color == 'fill'){
          gg = gg + scale_fill_brewer(palette = 'Set1')
        }
      }
      
    } else {
      if (n_batches <= 12){
        if(fill_or_color == 'color'){
          gg = gg + scale_color_brewer(palette = 'Set3')
        } else {
          if(fill_or_color == 'fill'){
            gg = gg + scale_fill_brewer(palette = 'Set3')
          }
        }
        
      } else {
        warning(sprintf('brewer palettes have maximally 12 colors, 
                        %s batches are specified,
                        consider defining color scheme with 
                        sample_annotation_to_colors function', 
                        n_batches))
      }
    }
    } else {
      #color vector provided by "sample_annotation_to_colors"
      if(fill_or_color == 'color'){
        gg = gg + scale_color_manual(values = color_scheme)
      } else {
        if(fill_or_color == 'fill'){
          gg = gg + scale_fill_manual(values = color_scheme)
        }
      }
    }
  if(fill_or_color == 'color'){
    gg = gg + labs(color = batch_col)
  } else {
    if(fill_or_color == 'fill'){
      gg = gg + labs(fill = batch_col)
    }
  }
  
  return(gg)
}

color_continuous <- function(color_scheme, batch_col, n_batches, fill_or_color, gg) {
  batch_vector = gg$data[[batch_col]]
  lab_datetime <- pretty(batch_vector)
  
  gg = gg + aes(color = as.numeric(!!sym(batch_col)))
  #Define the color scheme on the fly
  if(length(color_scheme) == 1 && color_scheme == 'brewer'){
    
    if(fill_or_color == 'color'){
      gg = gg + scale_color_distiller(palette = 'PiYG',
                                      breaks = as.numeric(lab_datetime), 
                                      labels = lab_datetime)+
        labs(color=batch_col)
    } else {
      if(fill_or_color == 'fill'){
        gg = gg + scale_color_distiller(palette = 'PiYG',
                                        breaks = as.numeric(lab_datetime), 
                                        labels = lab_datetime)+
          labs(fill=batch_col)
      }
    }
  } else {
    #color vector provided by "sample_annotation_to_colors"
    if(fill_or_color == 'color'){
      gg = gg + scale_color_gradientn(colors = color_scheme,
                                      breaks = as.numeric(lab_datetime), 
                                      labels = lab_datetime)+
        labs(color=batch_col)
    } else {
      if(fill_or_color == 'fill'){
        gg = gg + scale_color_gradientn(colors = color_scheme,
                                        breaks = as.numeric(lab_datetime), 
                                        labels = lab_datetime)+
          labs(fill=batch_col)
      }
    }
  }
  return(gg)
}

color_by_factor <- function(color_by_batch, batch_col, gg, color_scheme, 
                            sample_annotation, fill_or_color = 'color') {
  if(color_by_batch & !is.null(batch_col)){
    
    batch_vector <- sample_annotation[[batch_col]]
    n_batches <- length(unique(batch_vector))
    
    is_factor = is_batch_factor(batch_vector, color_scheme)
    
    is_numeric = (!is_factor) && 
      (is.numeric(batch_vector) || is.POSIXct(batch_vector))
    
    if (is_numeric && 
        (n_batches <= 10 || n_batches < 0.1*nrow(sample_annotation))){
      warning(sprintf('%s column has very few values, but is numeric-like,
                      should it be treated as factor? 
                      \nThen modify it with as.factor() function', batch_col))
    }
    
    if(is.null(color_scheme)){
      color_scheme = 'brewer'
    }
    
    if(is_factor){
      gg = color_discrete(color_scheme, batch_col, n_batches, fill_or_color, gg)
    } else if(is_numeric) {
      gg = color_continuous(color_scheme, batch_col, n_batches, fill_or_color, gg)
    } else{
      stop('batch_col class is neither factor-like, nor numeric-like,  
           check sample_annotation and/or color_scheme')
    }
    
    
    
    } else {
      if (is.null(batch_col)){
        stop('Coloring column not defined, please define the color column!')
      }
  }
  return(gg)
}
