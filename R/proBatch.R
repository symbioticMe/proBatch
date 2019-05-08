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
#' \code{browseVignettes(package = "proBatch")}
#'
#' @section Section:
#' Common arguments to the functions.
#' 
#' @param df_long data frame where each row is a single feature in a single
#'   sample. It minimally has a \code{sample_id_col}, a \code{feature_id_col} and a
#'   \code{measure_col}, but usually also an \code{m_score} (in OpenSWATH output result
#'   file)
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. Usually the log
#'   transformed version of the original data
#' @param sample_annotation data frame with: \enumerate{ \item \code{sample_id_col}
#'   (this can be repeated as row names) \item biological covariates \item
#'   technical covariates (batches etc) }
#' @param sample_id_col name of the column in \code{sample_annotation} table, where the
#'   filenames (colnames of the data matrix are found).
#' @param measure_col if \code{df_long} is among the parameters, it is the
#'   column with expression/abundance/intensity; otherwise, it is used
#'   internally for consistency.
#' @param feature_id_col name of the column with feature/gene/peptide/protein
#'   ID used in the long format representation \code{df_long}. In the wide
#'   formatted representation \code{data_matrix} this corresponds to the row
#'   names.
#' @param batch_col column in \code{sample_annotation} that should be used for
#'   batch comparison
#' @param order_col column in \code{sample_annotation} that determines sample order. It is
#'    used for in initial assessment plots (\link{plot_sample_mean_or_boxplot}) and 
#'    feature-level diagnostics (\link{feature_level_diagnostics}). Can be `NULL` 
#'    if sample order is irrelevant (e.g. in genomic experiments). For more details,
#'    order definition/inference, see \link{define_sample_order} and \link{date_to_sample_order}
#' @param facet_col column  in \code{sample_annotation} with a batch factor to separate 
#' plots into facets; usually 2nd to \code{batch_col}. Most meaningful for multi-instrument 
#' MS experiments (where each instrument has its own order-associated effects, see \code{order_col}) 
#' or simultaneous examination of two batch factors (e.g. preparation day and measurement day).
#' For single-instrument case should be set to `NULL`
#' @param color_by_batch (logical) whether to color points and connecting lines 
#' by batch factor as defined by \code{batch_col}.
#' @param qual_col column to color point by certain value denoted 
#' by \code{color_by_qual_value}. Design with inferred/requant values in openSWATH output data, 
#' which means argument value has to be set to \code{m_score}.
#' @param qual_value value in \code{qual_col} to color. For OpenSWATH data,
#' this argument value has to be set to \code{2} (this is an \code{m_score} value for requants).
#' @param plot_title title of the plot (usually, processing step + representation
#'   level (fragments, transitions, proteins))
#' @param theme ggplot theme, by default \code{classic}. Can be easily overriden 
#'
#' @import dplyr
#' @import ggfortify
#' @import ggplot2
#' @import reshape2
#' @import tibble
#' @import lazyeval
#' @importFrom corrplot corrplot.mixed
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices png pdf dev.off
#' @importFrom magrittr %>%
#' @importFrom pheatmap pheatmap
#' @importFrom purrr pmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom rlang :=
#' @importFrom rlang !!
#' @importFrom rlang !!!
#' @importFrom rlang sym syms
#' @importFrom tidyr complete nest unnest
#' @importFrom utils combn 
#' @importFrom scales brewer_pal
#' @importFrom stats as.formula complete.cases dist hclust loess median 
#' @importFrom stats model.matrix prcomp predict reformulate setNames
#' @importFrom viridis viridis_pal
#' @importFrom WGCNA plotDendroAndColors standardColors
#' 
#' @docType package
#' @name proBatch
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", "data",
  "batch_size", "batch_the_same", "batch_total", "category", "dateTime", 
  "fit", "label", "mean_fit", "median_batch", "median_global", 
  "median_run", "optimise_bw", "optimise_df", "peptide_col_name", 
  "sample_annotatation_col", "Step", "tipping.poings", "Var1", "Var2"))
NULL


