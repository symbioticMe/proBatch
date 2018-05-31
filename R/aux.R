#' Long to wide conversion
#'
#' Convert from a long data frame representation to a wide matrix representation
#'
#' @inheritParams proBatch
#'
#' @return \code{data_matrix} (\link{proBatch}) like matrix (features in rows, samples in columns)
#'
<<<<<<< HEAD
#' @export
#'
#' @family matrix manipulation functions
#'
=======
#' @family matrix manipulation functions
#'
#' @export
#' @import tibble
#' @importFrom magrittr %>%
#' @import dplyr
#' @import reshape2
#'
>>>>>>> f81d55e1139081fef18f4e9306967ce069471b96
convert_to_matrix <- function(df_long,
                              feature_id_column = 'peptide_group_label',
                              measure_column = 'Intensity',
                              sample_id_col = 'FullRunName') {
  casting_formula =  as.formula(paste(feature_id_column, sample_id_col,
                                      sep =  " ~ "))
  proteome_wide = dcast(df_long, formula = casting_formula,
                        value.var = measure_column) %>%
    column_to_rownames(feature_id_column) %>%
    as.matrix()
  return(proteome_wide)
}


#' Wide to long conversion
#'
#' Convert from wide matrix to a long data frame representation
#'
#' @inheritParams proBatch
#'
#' @param step normalization step (e.g. `Raw` or `Quantile_normalized` or
#'   `qNorm_ComBat`). Useful if consecutive steps are compared in plots. Note
#'   that in plots these are usually ordered alphabetically, so it's worth
#'   naming with numbers, e.g. `1_raw`, `2_quantile`
#'
#' @return \code{df_long} (\link{proBatch}) like data frame
#'
#' @family matrix manipulation functions
#'
#' @export
#'
matrix_to_long <- function(data_matrix, sample_annotation,
                           feature_id_column = 'peptide_group_label',
                           measure_column = 'Intensity',
                           sample_id_col = 'FullRunName',
                           step = NA){
  df_long = data_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = feature_id_column) %>%
    melt(id.var = feature_id_column, value.name = measure_column,
         variable.name = sample_id_col, factorsAsStrings = F) %>%
    mutate(Step = step) %>%
    merge(sample_annotation)
  return(df_long)
}


#' Join data matrices
#'
#' Joins 2 or more data matrices
#'
#' @inheritParams proBatch
#' @inheritParams matrix_to_long
#' @param matrices_list list of matrices in \code{data_matrix} (\link{proBatch}) format to be joined
#"
#' @return \code{df_long} (\link{proBatch}) like data frame with a row, having entries for:
#' \enumerate{
#'   \item \code{feature_id_col} (e.g. peptide name)
#'   \item \code{sample_id_col} (e.g. filename)
#'   \item \code{measure_col} (e.g. intensity/expression)
#'   \item \code{step} (e.g. `raw`, `quantile_norm`)
#' }
#'
#' @family matrix manipulation functions
#'
#' @export
#'
join_data_matrices <- function(matrices_list, step,
                               sample_annotation, measure_column = 'Intensity'){
  long_df_list = lapply(1:length(matrices_list), function(i){
    matrix_to_long(matrices_list[[i]], sample_annotation = sample_annotation,
                   measure_column = measure_column, step = step[i])
  })
  joined_df = do.call(rbind, long_df_list)

}
