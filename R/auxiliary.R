#' Long to wide conversion
#'
#' Convert from a long data frame representation to a wide matrix representation
#'
#' @inheritParams proBatch
#'
#' @return \code{data_matrix} (\link{proBatch}) like matrix (features in rows, samples in columns)
#'
#' @family matrix manipulation functions
#'
#' @export
#'
long_to_matrix <- function(df_long,
                           feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity',
                           sample_id_col = 'FullRunName') {
  casting_formula =  as.formula(paste(feature_id_col, sample_id_col,
                                      sep =  " ~ "))
  proteome_wide = dcast(df_long, formula = casting_formula,
                        value.var = measure_col) %>%
    column_to_rownames(feature_id_col) %>%
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
matrix_to_long <- function(data_matrix, sample_annotation = NULL,
                           feature_id_col = 'peptide_group_label',
                           measure_col = 'Intensity',
                           sample_id_col = 'FullRunName',
                           step = NULL){
  
  if(setequal(unique(sample_annotation[[sample_id_col]]), unique(colnames(data_matrix))) == FALSE){
    warning('Sample IDs in sample annotation not consistent with samples in input data.')}
  
  df_long = data_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = feature_id_col) %>%
    melt(id.var = feature_id_col, value.name = measure_col,
         variable.name = sample_id_col, factorsAsStrings = F)
  if(!is.null(step)){
    df_long = df_long %>%
      mutate(Step = step)
  }
  if(!is.null(sample_annotation))
    df_long = df_long %>%
      merge(sample_annotation, by = sample_id_col)
  return(df_long)
}



#' Create light-weight peptide annotation data frame for selection of illustrative proteins
#'
#' @param df_long
#' @param peptide_col column containing peptide ID
#' @param protein_col one or more columns contatining protein ID
#'
#' @return
#' @export
#'
#' @examples \donotrun{
#' peptide_annotation = create_peptide_annotation(example_proteome)
#' peptide_summary =
#' }
#'
#' @seealso \code{\link{plot_peptides_of_one_protein}}, \code{\link{plot_corrplot_protein}},
#' \code{\link{plot_within_prot_distribution}}
create_peptide_annotation <- function(df_long, feature_id_col = 'peptide_group_label',
                                      annotation_col = c('RT', 'Intensity', "ProteinName" )){
  peptide_annotation = df_long %>%
    select(one_of(c(feature_id_col, annotation_col))) %>%
    distinct()
  return(peptide_annotation)
}
