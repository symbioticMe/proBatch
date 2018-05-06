#' Auxiliary functions for data matrix manipulation
#'
#'  Converting from long to wide (matrix), wide (matrix) to long, joining the data matrices
#'
#' @param df_long data frame where each row is a single feature in a single sample,
#' thus it has minimally, `sample_id_col`, `feature_id_column` and `measure_column`, but usually also `m_score` (in OpenSWATH output result file)
#' @param data_matrix features (in rows) vs samples (in columns) matrix,
#' with feature IDs in rownames and file/sample names as colnames.
#' Usually the log transformed version of the original data
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names) 2) biological and 3) technical covariates (batches etc)
#' @param sample_id_col name of the column in sample_annotation file,
#' where the filenames (colnames of the data matrix are found)
#' @param measure_column if `df_long` is among the parameters, it is the column with expression/abundance/intensity,
#' otherwise, it is used internally for consistency
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be repeated as row names) 2) biological and 3) technical covariates (batches etc)
#' @param feature_id_column name of the column with feature/gene/peptide/protein ID used with long format matrix (`df_long`). In wide format (`data_matrix`) this would be the row name
#' @param step normalization step (e.g. `Raw` or `Quantile_normalized` or `qNorm_ComBat`).
#' Useful if consecutive steps are compared in plots.
#' Note that in plots these are usually ordered alphabetically, so it's worth naming with numbers, e.g. `1_raw`, `2_quantile`
#' @name aux_config

#' @name aux_config
#'
#' @return `data_matrix`-like matrix (features in rows, samples in columns)
#' @export
#' @import tibble
#' @importFrom magrittr %>%
#' @import dplyr
#' @import reshape2
#'
#' @examples
convert_to_matrix <- function(df_long,
                              feature_id_column = 'peptide_group_label',
                              measure_column = 'Intensity',
                              sample_id_col = 'FullRunName'){
  casting_formula =  as.formula(paste(feature_id_column, sample_id_col,
                                      sep =  " ~ "))
  proteome_wide = dcast(df_long, formula=casting_formula,
                        value.var=measure_column) %>%
    column_to_rownames(feature_id_column) %>%
    as.matrix()
  return(proteome_wide)
}

#' @name aux_config
#'
#' @return `df_long`-like object
#' @export
#' @importFrom magrittr %>%
#' @import dplyr
#' @import reshape2
#'
#' @examples
matrix_to_long <- function(data_matrix, sample_annotation,
                           feature_id_column = 'peptide_group_label',
                           measure_col = 'Intensity', step = NA){
  df_long = data_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = feature_id_column) %>%
    melt(id.var = 'peptide_group_label', value.name = measure_col,
         variable.name = 'FullRunName', factorsAsStrings = F) %>%
    mutate(Step = step) %>%
    merge(sample_annotation)
  return(df_long)
}

#' @param matrix_list list of matrices in `data_matrix` format
#' @name aux_config
#'
#' @return `df_long` like data frame with a row, having entries for
#' 1) `feature_id_col` (e.g. peptide name), 2) `sample_id_col` (e.g. filename),
#' 3) `measure_col` (e.g. intensity/expression) and 4) `step` (e.g. `raw`, `quantile_norm`)
#' @export
#'
#' @examples
join_data_matrices <- function(matrix_list, step,
                               sample_annotation, measure_column = 'Intensity'){
  long_df_list = lapply(1:length(matrix_list), function(i){
    matrix_to_long(matrix_list[[i]], sample_annotation = sample_annotation,
                   measure_col = measure_column, step = step[i])
  })
  joined_df = do.call(rbind, long_df_list)

}
