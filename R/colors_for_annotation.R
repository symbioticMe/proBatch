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
  
  number_colors_for_factors = sum(nlev_covariate)
  standard_colors_base = grep('(white|(gr(a|e)y))', standardColors(45), 
                              value = TRUE, invert = TRUE)
  
  qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
  brewer_colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                                rownames(qual_col_pals)))
  if(number_colors_for_factors <= 45){
    colors = standard_colors_base
  } else 
    if(number_colors_for_factors <= length(brewer_colors)){
      colors = brewer_colors
    } else if (number_colors_for_factors <= sum(length(brewer_colors), 45)){
      colors = c(standard_colors_base, brewer_colors)
    } else {
      colors_to_sample = grep('(white|(gr(a|e)y))', standardColors(), 
                              value = TRUE, invert = TRUE)
      #equivalent to set.seed(1); indices_random = sample(1:length(colors_to_sample))
      indices_random = c(270L, 44L, 381L, 308L, 271L, 263L, 167L, 215L, 377L, 349L, 
                         26L, 141L, 19L, 57L, 248L, 137L, 240L, 406L, 107L, 252L, 20L, 
                         145L, 75L, 176L, 15L, 347L, 405L, 204L, 364L, 131L, 99L, 231L, 
                         352L, 331L, 201L, 357L, 49L, 192L, 18L, 14L, 402L, 188L, 211L, 
                         298L, 302L, 150L, 265L, 338L, 59L, 130L, 239L, 309L, 33L, 114L, 
                         129L, 353L, 284L, 264L, 103L, 135L, 354L, 190L, 158L, 175L, 379L, 
                         395L, 180L, 243L, 210L, 36L, 2L, 123L, 203L, 209L, 83L, 340L, 
                         157L, 25L, 407L, 189L, 350L, 177L, 132L, 51L, 320L, 368L, 341L, 
                         101L, 170L, 411L, 316L, 296L, 40L, 45L, 333L, 366L, 251L, 376L, 
                         179L, 152L, 115L, 56L, 138L, 134L, 370L, 392L, 288L, 393L, 322L, 
                         116L, 286L, 156L, 326L, 165L, 276L, 151L, 186L, 299L, 117L, 306L, 
                         382L, 219L, 220L, 93L, 267L, 412L, 312L, 147L, 244L, 184L, 386L, 
                         212L, 30L, 408L, 260L, 222L, 360L, 314L, 62L, 327L, 4L, 60L, 
                         24L, 345L, 301L, 378L, 330L, 323L, 78L, 185L, 67L, 257L, 98L, 
                         329L, 64L, 124L, 390L, 344L, 69L, 153L, 187L, 361L, 358L, 398L, 
                         82L, 221L, 321L, 399L, 295L, 66L, 171L, 85L, 54L, 198L, 79L, 
                         196L, 237L, 249L, 52L, 242L, 258L, 58L, 290L, 31L, 415L, 245L, 
                         317L, 403L, 63L, 154L, 266L, 38L, 74L, 384L, 9L, 401L, 193L, 
                         268L, 13L, 283L, 206L, 334L, 388L, 22L, 277L, 88L, 272L, 385L, 
                         339L, 37L, 343L, 168L, 346L, 235L, 126L, 241L, 108L, 112L, 225L, 
                         16L, 373L, 89L, 182L, 118L, 304L, 213L, 214L, 410L, 208L, 166L, 
                         351L, 285L, 84L, 92L, 128L, 169L, 297L, 292L, 5L, 146L, 71L, 
                         375L, 139L, 133L, 200L, 261L, 418L, 10L, 127L, 289L, 259L, 122L, 
                         233L, 3L, 53L, 383L, 318L, 104L, 396L, 161L, 155L, 371L, 61L, 
                         282L, 7L, 149L, 11L, 254L, 199L, 73L, 232L, 72L, 76L, 342L, 1L, 
                         355L, 105L, 397L, 12L, 32L, 90L, 262L, 315L, 77L, 335L, 325L, 
                         313L, 207L, 121L, 374L, 110L, 348L, 46L, 394L, 159L, 367L, 238L, 
                         380L, 229L, 365L, 48L, 274L, 387L, 311L, 273L, 409L, 23L, 413L, 
                         217L, 226L, 303L, 109L, 162L, 224L, 172L, 356L, 250L, 140L, 174L, 
                         86L, 400L, 269L, 143L, 234L, 337L, 8L, 195L, 281L, 113L, 29L, 
                         300L, 389L, 163L, 255L, 362L, 47L, 65L, 287L, 216L, 80L, 246L, 
                         28L, 278L, 100L, 43L, 191L, 87L, 35L, 391L, 404L, 164L, 68L, 
                         136L, 102L, 119L, 50L, 328L, 81L, 205L, 332L, 144L, 120L, 372L, 
                         310L, 307L, 39L, 280L, 230L, 173L, 227L, 106L, 96L, 42L, 178L, 
                         34L, 202L, 336L, 236L, 194L, 275L, 21L, 253L, 294L, 91L, 369L, 
                         223L, 148L, 319L, 416L, 419L, 293L, 94L, 183L, 70L, 6L, 111L, 
                         279L, 324L, 197L, 160L, 17L, 218L, 142L, 256L, 247L, 228L, 181L, 
                         27L, 291L, 55L, 41L, 359L, 95L, 417L, 125L, 363L, 305L, 414L, 
                         97L)
      standard_colors_all = colors_to_sample[indices_random]
      colors = standard_colors_all
  }
  
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
#' @param column column of the data whose rare categories need to be merged to 
#' "other"
#' @param rare_thr minimal number of times for a category to be represented to 
#' be declared as "rare" and converted to "other"
#'
#' @keywords internal
#'
#' @return column with rare occurrences replaced by other
#' 
#' @examples 
#' column <- factor(c("A", "B", "A", "C", "D", "D", "E"))
#' merge_rare_levels(column, rare_thr = 2)
#' # [1] A     other A     other D     D     other
#' # Levels: A D other
#' @export
merge_rare_levels <- function(column, rare_thr = 2) {
  is_factor_col = is.factor(column)
  tb_col = table(column)
  if (is_factor_col)
    column = as.character(column)
  column[column %in% names(tb_col)[tb_col < rare_thr]] = 'other'
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
                                        numeric_columns = c('DateTime','order'),
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
    warning('numeric columns not specified, 
            extracting numeric columns from factors')
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
    is_not_factor <- vapply(factor_columns, 
                            function(col) is.numeric(sample_annotation[[col]]),
                            FUN.VALUE = logical(1))
    numeric_columns = factor_columns[is_not_factor]
    numeric_columns = c(numeric_columns, date_columns)
    for (numcol in numeric_columns){
      batch_vector = sample_annotation[[numcol]]
      n_batches = length(unique(batch_vector))
      if (n_batches <= 10 || n_batches < 0.1*nrow(sample_annotation)){
        warning(sprintf('%s column has very few values, but is numeric-like,
                        should it be treated as factor? 
                        \nuse both factor_columns and 
                        numeric_columns parameters', numcol))
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
        mutate_if(check_rare_levels, list(~ merge_rare_levels(.)))
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


color_discrete <- function(color_scheme, batch_col, n_batches, fill_or_color, 
                           gg) {
  
  if(fill_or_color == 'color'){
    gg = gg + aes(color = as.factor(!!sym(batch_col)))
  } else {
    if(fill_or_color == 'fill'){
      gg = gg + aes(fill = as.factor(!!sym(batch_col)))
    }
  }
  
  #Define the color scheme on the fly
  gg = add_color_scheme_discrete(color_scheme, n_batches, fill_or_color, 
                                 gg, batch_col)
  
  return(gg)
}

color_continuous <- function(color_scheme, batch_col, n_batches, 
                             fill_or_color, gg) {
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
    
    if((length(color_scheme) == 1) && color_scheme == 'brewer'){
      warning('color_scheme will be inferred automatically.
              Numeric/factor columns are guessed, for more controlled color 
              mapping use sample_annotation_to_colors()')
    }
    
    if(is_factor){
      gg = color_discrete(color_scheme, batch_col, n_batches, fill_or_color, gg)
    } else if(is_numeric) {
      gg = color_continuous(color_scheme, batch_col, n_batches, fill_or_color, 
                            gg)
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
