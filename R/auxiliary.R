#' Long to wide data format conversion
#'
#' Convert from a long data frame representation to a wide matrix representation
#'
#' @inheritParams proBatch
#'
#' @return \code{data_matrix} (\link{proBatch}) like matrix
#' (features in rows, samples in columns)
#'
#' @family matrix manipulation functions
#' @examples
#' data("example_proteome", package = "proBatch")
#' proteome_matrix <- long_to_matrix(example_proteome)
#'
#' @export
#'
long_to_matrix <- function(df_long,
                           feature_id_col = "peptide_group_label",
                           measure_col = "Intensity",
                           sample_id_col = "FullRunName",
                           qual_col = NULL,
                           qual_value = 2) {
  casting_formula <- as.formula(paste(
    feature_id_col, sample_id_col,
    sep = " ~ "
  ))
  if (!is.null(qual_col)) {
    message("removing imputed values (requants)")
    df_long <- df_long %>%
      mutate(!!sym(measure_col) := ifelse(
        !!sym(qual_col) == qual_value, NA, !!sym(measure_col)
      ))
  }
  proteome_wide <- dcast(
    df_long,
    formula = casting_formula,
    value.var = measure_col
  ) %>%
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
#' @param step normalization step (e.g. \code{Raw} or \code{Normalized}.
#' Useful if consecutive steps are compared in plots. Note
#'   that in plots these are usually ordered alphabetically, so it's worth
#'   naming with numbers, e.g. \code{1_raw}, \code{2_quantile}
#'
#' @return \code{df_long} (\link{proBatch}) like data frame
#'
#' @family matrix manipulation functions
#' @examples
#' # Load necessary datasets
#' data(
#'   list = c("example_sample_annotation", "example_proteome_matrix"),
#'   package = "proBatch"
#' )
#' # Convert matrix to long format
#' proteome_long <- matrix_to_long(
#'   example_proteome_matrix,
#'   example_sample_annotation
#' )
#'
#' @export
#'
matrix_to_long <- function(data_matrix, sample_annotation = NULL,
                           feature_id_col = "peptide_group_label",
                           measure_col = "Intensity",
                           sample_id_col = "FullRunName",
                           step = NULL) {
  df_long <- data_matrix %>%
    as.data.frame() %>%
    rownames_to_column(var = feature_id_col) %>%
    melt(
      id.var = feature_id_col, value.name = measure_col,
      variable.name = sample_id_col, factorsAsStrings = FALSE
    )
  if (!is.null(step)) {
    df_long <- df_long %>%
      mutate(Step = step)
  }

  if (!is.null(sample_annotation)) {
    df_long <- check_sample_consistency(
      sample_annotation = sample_annotation,
      sample_id_col = sample_id_col,
      df_long = df_long,
      batch_col = NULL, order_col = NULL,
      facet_col = NULL, merge = FALSE
    )
  }

  return(df_long)
}



#' Prepare peptide annotation from long format data frame
#'
#' Create light-weight peptide annotation data frame
#' for selection of illustrative proteins
#'
#' @inheritParams proBatch
#'
#' @return data frame containing petpide annotations
#' @export
#' @examples
#' data("example_proteome", package = "proBatch")
#' generated_peptide_annotation <- create_peptide_annotation(
#'   example_proteome,
#'   feature_id_col = "peptide_group_label",
#'   protein_col = c("Protein")
#' )
#'
#' @seealso \code{\link{plot_peptides_of_one_protein}},
#' \code{\link{plot_protein_corrplot}}
create_peptide_annotation <- function(df_long,
                                      feature_id_col = "peptide_group_label",
                                      protein_col = c("ProteinName", "Gene")) {
  if (!all(protein_col %in% names(df_long))) {
    stop(
      sprintf("Column %s is not in the data"),
      setdiff(names(df_long), protein_col)
    )
  }
  peptide_annotation <- df_long %>%
    select(one_of(c(feature_id_col, protein_col))) %>%
    distinct()
  return(peptide_annotation)
}
