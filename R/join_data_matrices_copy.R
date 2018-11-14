
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
join_data_matrices <- function(matrices_list, step = NULL,
                               sample_annotation = NULL,
                               feature_id_col = 'peptide_group_label',
                               measure_col = 'Intensity',
                               sample_id_col = 'FullRunName'){
  if(is.null(step)) step = names(matrices_list)
  if(length(step) != length(matrices_list)){
    stop('check the step names and list of matrices')
  }
  long_df_list = lapply(1:length(matrices_list), function(i){
    matrix_to_long(matrices_list[[i]], sample_annotation = sample_annotation,
                   measure_col = measure_col, step = step[i],
                   feature_id_col = feature_id_col, sample_id_col = sample_id_col)
  })
  joined_df = do.call(rbind, long_df_list)
  return(joined_df)
}