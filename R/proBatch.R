#' proBatch: A package for diagnostics and correction of batch effects,
#' primarily in proteomics
#'
#' The proBatch package contains functions for analyzing and correcting batch effects 
#' and other unwanted technical variation from high-thoughput experiments. Although 
#' the package has primarily been developed for mass spectrometry proteomics (DIA/SWATH),
#' it should also be applicable to most omic data with minor adaptations.
#' It addresses the following needs: \itemize{ \item prepare the data for analysis
#'\item Visualize batch effects in sample-wide and feature-level;
#' \item Normalize and correct for batch effects.
#' }
#'
#' To learn more about proBatch, start with the vignettes:
#' `browseVignettes(package = "proBatch")`
#'
#' @section Section:
#' Common arguments to the functions.
#' 
#'
#' @param df_long data frame where each row is a single feature in a single
#'   sample. It minimally has a \code{sample_id_col}, a \code{feature_id_col} and a
#'   \code{measure_col}, but usually also an \code{m_score} (in OpenSWATH output result
#'   file)
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data matrix with: \enumerate{ \item \code{sample_id_col}
#'   (this can be repeated as row names) \item biological covariates \item
#'   technical covariates (batches etc) }
#' @param sample_id_col name of the column in sample_annotation file, where the
#'   filenames (colnames of the data matrix are found)
#' @param batch_col column in \code{sample_annotation} that should be used for
#'   batch comparison
#' @param order_col column in \code{sample_annotation} that determines sample order. It is
#'    used for certain diagnostics and normalisations.
#' @param measure_col if \code{df_long} is among the parameters, it is the
#'   column with expression/abundance/intensity; otherwise, it is used
#'   internally for consistency
#' @param feature_id_col name of the column with feature/gene/peptide/protein
#'   ID used in the long format representation \code{df_long}. In the wide
#'   formatted representation \code{data_matrix} this corresponds to the row
#'   names.
#' @param plot_title Title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param theme ggplot theme, by default `classic`. Can be easily overriden 
#'
#' @import dplyr
#' @import ggfortify
#' @import ggplot2
#' @import pheatmap
#' @import reshape2
#' @import tibble
#' @import lazyeval
#' @import WGCNA
#' @import viridis
#' @import RColorBrewer
#' @import readr
#' @import tidyverse
#' @importFrom corrplot corrplot.mixed
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom rlang UQ sym syms
#' @importFrom tidyr complete nest unnest
#' @importFrom data.table setDT
#' @importFrom data.table transpose
#' @importFrom data.table :=
#' @importFrom data.table IDateTime
#' @importFrom purrr map2
#' 
#' @docType package
#' @name proBatch
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "batch_size", "batch_the_same", "batch_total", "category", "dateTime", "fit", "label", "mean_fit", "median_batch", "median_global", "median_run", "optimise_bw", "optimise_df", "peptide_col_name", "sample_annotatation_col", "Step", "tipping.poings", "Var1", "Var2"))
NULL


