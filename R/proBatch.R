#' proBatch: A package for diagnostics and correction of batch effects,
#' primarily in proteomics
#'
#' It adresses the following needs: \itemize{ \item prepare the original data
#' (e.g. OpenSWATH output matrix and sample annotation file) for analysis.
#' However, you might need to use `SWATH2stats` additionally \item Diagnose
#' batch effects, sample-wide and feature-level \item Correct for batch effects
#' (normalize the data). Other useful package for this purpose is `Normalyzer`.
#' }
#'
#' To learn more about proBatch, start with the vignettes:
#' `browseVignettes(package = "proBatch")`
#'
#' @section common arguments to the functions:
#'
#' @param df_long data frame where each row is a single feature in a single
#'   sample. It minimally has a \code{sample_id_col}, a \code{feature_id_column} and a
#'   \code{measure_column}, but usually also an \code{m_score} (in OpenSWATH output result
#'   file)
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data matrix with: \enumerate{ \item \code{sample_id_col}
#'   (this can be repeated as row names) \item biological covariates \item
#'   technical covariates (batches etc) }
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param batch_column column in \code{sample_annotation} that should be used for
#'   batch comparison
#' @param order_column column in \code{sample_annotation} that determines sample order. It is
#'    used for certain diagnostics and normalisations.
#' @param measure_column if \code{df_long} is among the parameters, it is the
#'   column with expression/abundance/intensity; otherwise, it is used
#'   internally for consistency
#' @param feature_id_column name of the column with feature/gene/peptide/protein
#'   ID used in the long format representation \code{df_long}. In the wide
#'   formatted representation \code{data_matrix} this corresponds to the row
#'   names.
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param theme ggplot theme, by default `classic`. Can be easily overriden (see
#'   examples)
#'
#' @import dplyr
#' @import ggfortify
#' @import ggplot2
#' @import pheatmap
#' @import reshape2
#' @import tibble
#' @importFrom corrplot corrplot.mixed
#' @importFrom magrittr %>%
#' @importFrom  purrr map
#' @importFrom rlang UQ sym syms
#' @importFrom tidyr complete nest unnest
#'
#' @docType package
#' @name proBatch
NULL


# TODO: refactor sample_id_col to sample_id_column for naming consistency

# TODO: Pick a variable naming stlye and refactor the rest out (at the moment we have a_b, aB, and a.b)
