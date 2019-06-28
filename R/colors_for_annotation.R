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
  
  if (any(nlev_covariate > 20)) {
    warning('Some colors will be hard to distinguish\n')
  }
  if (any(nlev_covariate > 50)) {
    warning('Too many colors, consider merging some covariate 
                values for better visualisation\n')
  }

  colors = grep('(white|(gr(a|e)y)', standardColors(), value = T, invert = T)
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

map_numbers_to_colors <- function(annotation_df_numbers,
                                  palette_type = 'brewer') {
  n_colors_to_create <- ncol(annotation_df_numbers)
  if ((n_colors_to_create > 4 & palette_type == 'viridis')) {
    warning('Too many colors for viridis palette, 
            switching to Brewer palettes')
  }
  if ((n_colors_to_create > 18)) {
    stop('Not enough color paletters to visualize the annotation')
  }
  
  color_list <- lapply(seq_len(ncol(annotation_df_numbers)),
                       function(i)
                         generate_colors_for_numeric(i = i,
                                                     palette_type = palette_type))
  names(color_list) = names(annotation_df_numbers)
  return(color_list)
}


#' Generates color vector from continous palette
#'
#' Generates a vector of colors for a vector of numeric, POSIXct (i.e. the
#' (signed) number of seconds since the beginning of 1970 , or factors
#'
#' @param palette_type 'brewer' or 'viridis'
#' @param i if \code{palette_type} is 'brewer' the palette argument to
#'   \code{brewer_pal}. If \code{palette_type} is 'viridis' the option argument
#'   to \code{viridis_pal}
#'
#' @return vector of colors
#' @keywords internal
#' 
generate_colors_for_numeric <- function(palette_type = 'brewer',
                                        i = 1) {
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
  
  return(color_for_column)
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
#' @param numeric_columns columns of \code{sample_annotation} to be 
#' treated as continuous numeric values. 
#' @param rare_categories_to_other if \code{True} rare categories 
#' will be merged into the value \code{"other"}
#' @param guess_factors whether attempt which of the \code{factor_columns}
#'  are actually numeric
#' @param numeric_palette_type palette to be used for 
#' numeric values coloring (can be \code{'brewer' and 'viridis'})
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
#' numeric_columns = c('DateTime', 'order'))
#' @export
#' 
#' @name sample_annotation_to_colors
sample_annotation_to_colors <- function(sample_annotation,
                                        sample_id_col = 'FullRunName',
                                        factor_columns = c('MS_batch','EarTag', 
                                                           'digestion_batch',
                                                           "Strain", "Diet"),
                                        numeric_columns = c('DateTime', 'order'),
                                        rare_categories_to_other = TRUE,
                                        guess_factors = FALSE,
                                        numeric_palette_type = 'brewer') {
  
  sample_annotation = as.data.frame(sample_annotation)
  
  columns_for_color_mapping = union(factor_columns, numeric_columns)
  undefined_cols <- setdiff(names(sample_annotation), 
                            c(columns_for_color_mapping, sample_id_col))
  if (length(undefined_cols) > 0){
    warning(paste(c('The following columns will not be mapped to colors:',
                    undefined_cols, 'if these have to be mapped, please assign 
                    them to factor, date or numeric'), collapse = ' '))
  }
  
  sample_annotation = sample_annotation %>%
    select(one_of(columns_for_color_mapping))
  
  if (!is.null(numeric_columns)) {
    factor_columns = setdiff(factor_columns, numeric_columns)
    column_intersection <- intersect(factor_columns, numeric_columns)
    if (length(column_intersection) > 0) {
      warning(paste(c('The following columns are repeatedly listed among factors
                      and numeric-like variables:', column_intersection, 
                      '; they will be excluded from factors and mapped to 
                      continuous palettes'), 
                    collapse = ' '))
    }
  }
  
  if(is.null(numeric_columns) && guess_factors){
    warning('numeric columns not specified, extracting numeric columns from factors')
    which_factors = vapply(factor_columns, function(col) {
      batch_vector = sample_annotation[[col]]
      is_factor = is_batch_factor(batch_vector, color_scheme = NULL)
      return(is_factor)
    }, FUN.VALUE = logical(1))
    factor_candidates = factor_columns[which_factors]
    
    which_dates = vapply(factor_columns, 
                         function(col) is.POSIXct(sample_annotation[[col]]), 
                         FUN.VALUE = logical(1))
    date_columns = factor_columns[which_dates]
    numeric_columns = factor_columns[vapply(factor_columns, 
                                            function(col) is.numeric(sample_annotation[[col]]), 
                                            FUN.VALUE = logical(1))]
    numeric_columns = c(numeric_columns, date_columns)
    for (numcol in numeric_columns){
      batch_vector = sample_annotation[[numcol]]
      n_batches = length(unique(batch_vector))
      if (n_batches <= 10 || n_batches < 0.1*nrow(sample_annotation)){
        warning(sprintf('%s column has very few values, but is numeric-like,
                        should it be treated as factor? 
                        \nuse both factor_columns and numeric_columns parameters', numcol))
      }
    }
    factor_columns = factor_candidates
  }

  message('converting columns to corresponding classes 
          (factor, numeric)')
  sample_annotation <- sample_annotation %>%
    mutate_at(vars(factor_columns), as.factor) %>%
    mutate_at(vars(numeric_columns), as.numeric)
  
  if (!is.null(factor_columns) || length(factor_columns) != 0) {
    factor_df = sample_annotation %>%
      select(one_of(factor_columns))
    if (rare_categories_to_other) {
      factor_df = factor_df %>%
        mutate_if(check_rare_levels, funs(merge_rare_levels(.)))
    }
    list_of_col_for_factors = map_factors_to_colors(factor_df)
  } else {
    list_of_col_for_factors = list()
  }
  
  non_factor_cols = setdiff(columns_for_color_mapping, factor_columns)
  
  if (!is.null(non_factor_cols) &
      !identical(non_factor_cols, character(0))) {
    numeric_df = sample_annotation %>%
      select(one_of(non_factor_cols))
    list_of_col_for_numeric = map_numbers_to_colors(
      numeric_df, palette_type = numeric_palette_type)
  } else {
    list_of_col_for_numeric = list()
  }
  
  color_list = c(list_of_col_for_factors, list_of_col_for_numeric)
  
  return(color_list)
}

map_numeric_colors_to_intervals <- function(color_vector, col_values){
  breaks = pretty(col_values)
  col_intervals =cut(col_values, breaks = breaks)
  col_for_colors = colorRampPalette(color_vector)(nlevels(col_intervals))
  names(col_for_colors) = levels(col_intervals)
  col_colors = col_for_colors[col_intervals]
  return(col_colors)
}

#' Color list to data frame
#'
#' Turn color list to df (to use in the hierarchical clustering)
#'
#' @param color_list list of colors
#' @param sample_annotation factor-based configuration 
#' of the sample annotation
#' 
#' @return a data frame representation of the input color list
#' 
#' @keywords internal
#'
color_list_to_df <- function(color_list, sample_annotation, 
                             sample_id_col = 'FullRunName') {
  factors_to_map = intersect(names(sample_annotation), names(color_list))
  if(!setequal(names(sample_annotation), names(color_list))){
    warning('color list and sample annotation have different factors, 
            using only intersection in color scheme!')
  }
  list_df = lapply(factors_to_map, function(col_name) {
    
    col_values = sample_annotation[[col_name]]
    color_scheme <- color_list[[col_name]]
    is_factor = is_batch_factor(col_values, color_scheme=color_scheme)
    if(is_factor){
      col_colors = color_scheme[col_values]
    } else{
      col_colors = map_numeric_colors_to_intervals(color_scheme, col_values)
    }
    if(any(is.na(col_values))){
      col_colors[is.na(col_values)] = 'white'
    }
    return(col_colors)
  })
  names(list_df) = factors_to_map
  color_df = as.data.frame(do.call(cbind, list_df))
  rownames(color_df) = sample_annotation[[sample_id_col]]
  return(color_df)
}

add_color_scheme_discrete <- function(color_scheme, n_batches, fill_or_color, 
                                      gg, batch_col) {
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


color_discrete <- function(color_scheme, batch_col, n_batches, fill_or_color, gg) {
  
  if(fill_or_color == 'color'){
    gg = gg + aes(color = as.factor(!!sym(batch_col)))
  } else {
    if(fill_or_color == 'fill'){
      gg = gg + aes(fill = as.factor(!!sym(batch_col)))
    }
  }
  
  #Define the color scheme on the fly
  gg = add_color_scheme_discrete(color_scheme, n_batches, fill_or_color, gg, batch_col)
  
  return(gg)
}

color_continuous <- function(color_scheme, batch_col, n_batches, fill_or_color, gg) {
  batch_vector = gg$data[[batch_col]]
  lab_datetime <- pretty(batch_vector)
  
  if(fill_or_color == 'color'){
    gg = gg + aes(color = as.numeric(!!sym(batch_col)))
  } else {
    if(fill_or_color == 'fill'){
      gg = gg + aes(fill = as.numeric(!!sym(batch_col)))
    }
  }
  
  #Define the color scheme on the fly
  if(length(color_scheme) == 1 && color_scheme == 'brewer'){
    
    if(fill_or_color == 'color'){
      gg = gg + scale_color_distiller(palette = 'PiYG',
                                      breaks = as.numeric(lab_datetime), 
                                      labels = lab_datetime)+
        labs(color=batch_col)
    } else {
      if(fill_or_color == 'fill'){
        gg = gg + scale_fill_distiller(palette = 'PiYG',
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
        gg = gg + scale_fill_gradientn(colors = color_scheme,
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
    
    if(color_scheme == 'brewer'){
      warning('color_scheme will be inferred automatically.
              Numeric/factor columns are guessed, for more controlled color mapping use 
              sample_annotation_to_colors()')
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
