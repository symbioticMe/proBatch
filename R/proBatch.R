#' proBatch: A package for diagnostics and correction of batch effects,
#' primarily in proteomics
#'
#' The proBatch package contains functions for analyzing and correcting batch 
#' effects (unwanted technical variation) from high-thoughput experiments. 
#' Although the package has primarily been developed for mass spectrometry 
#' proteomics (DIA/SWATH), it has been designed be applicable to most omic data 
#' with minor adaptations.
#' It addresses the following needs: 
#' \itemize{ \item prepare the data for analysis
#' \item Visualize batch effects in sample-wide and feature-level;
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
#'   sample. It minimally has a \code{sample_id_col}, a \code{feature_id_col} 
#'   and a \code{measure_col}, but usually also an \code{m_score} (in OpenSWATH 
#'   output result file). See \code{help("example_proteome")} for more details.
#' @param data_matrix features (in rows) vs samples (in columns) matrix, with
#'   feature IDs in rownames and file/sample names as colnames. 
#'   See "example_proteome_matrix" for more details (to call the description, 
#'   use \code{help("example_proteome_matrix")})
#' @param sample_annotation data frame with: 
#' \enumerate{ \item \code{sample_id_col} (this can be repeated as row names) 
#'   \item biological covariates 
#'   \item technical covariates (batches etc) }. 
#'   See \code{help("example_sample_annotation")}
#' @param sample_id_col name of the column in \code{sample_annotation} table, 
#' where the filenames (colnames of the \code{data_matrix} are found).
#' @param measure_col if \code{df_long} is among the parameters, it is the
#'   column with expression/abundance/intensity; otherwise, it is used
#'   internally for consistency.
#' @param feature_id_col name of the column with feature/gene/peptide/protein
#'   ID used in the long format representation \code{df_long}. In the wide
#'   formatted representation \code{data_matrix} this corresponds to the row
#'   names.
#' @param batch_col column in \code{sample_annotation} that should be used for
#' batch comparison (or other, non-batch factor to be mapped to color in plots).
#' @param order_col column in \code{sample_annotation} that determines sample 
#' order. It is used for in initial assessment plots 
#' (\link{plot_sample_mean_or_boxplot}) and  feature-level diagnostics 
#' (\link{feature_level_diagnostics}). Can be `NULL` 
#'    if sample order is irrelevant (e.g. in genomic experiments). For more 
#'    details,
#'    order definition/inference, see \link{define_sample_order} and 
#'    \link{date_to_sample_order}
#' @param facet_col column  in \code{sample_annotation} with a batch factor to 
#' separate plots into facets; usually 2nd to \code{batch_col}. Most meaningful 
#' for multi-instrument MS experiments (where each instrument has its own 
#' order-associated effects (see \code{order_col}) or simultaneous examination 
#' of two batch factors (e.g. preparation day and measurement day).
#' For single-instrument case should be set to `NULL`
#' @param color_by_batch (logical) whether to color points and connecting lines 
#' by batch factor as defined by \code{batch_col}.
#' @param peptide_annotation long format data frame with peptide ID and their 
#' corresponding protein and/or gene annotations. 
#' See \code{help("example_peptide_annotation")}.
#' @param color_scheme a named vector of colors to map to \code{batch_col}, 
#' names corresponding to the levels of the factor. For continuous variables, 
#' vector doesn't need to be named.
#' @param color_list list, as returned by \code{sample_annotation_to_colors}, 
#' where each item contains a color vector for each factor to be mapped to the 
#' color.
#' @param factors_to_plot vector of technical and biological covariates to be 
#' plotted in this diagnostic plot (assumed to be present in 
#' \code{sample_annotation})
#' @param protein_col column where protein names are specified
#' @param no_fit_imputed (logical) whether to use imputed (requant) values, as flagged in 
#' \code{qual_col} by \code{qual_value} for data transformation
#' @param qual_col column to color point by certain value denoted 
#' by \code{qual_value}. Design with inferred/requant values in 
#' OpenSWATH output data, 
#' which means argument value has to be set to \code{m_score}.
#' @param qual_value value in \code{qual_col} to color. For OpenSWATH data,
#' this argument value has to be set to \code{2} (this is an \code{m_score} 
#' value for imputed values (requant values).
#' @param plot_title title of the plot (e.g., processing step + representation
#'   level (fragments, transitions, proteins) + purpose (meanplot/corrplot etc))
#' @param keep_all when transforming the data (normalize, correct) - acceptable 
#' values: all/default/minimal (which set of columns be kept).
#' @param theme ggplot theme, by default \code{classic}. Can be easily overriden
#' @param filename path where the results are saved. 
#' If null the object is returned to the active window;
#' otherwise, the object is save into the file. Currently only 
#' pdf and png format is supported
#' @param width option  determining the output image width
#' @param height option  determining the output image width
#' @param units units: 'cm', 'in' or 'mm'
#' @param base_size base font size
#'
#' @import dplyr
#' @import ggfortify
#' @import ggplot2
#' @import reshape2
#' @import lazyeval
#' @importFrom corrplot corrplot.mixed
#' @importFrom grDevices colorRampPalette
#' @importFrom grDevices png pdf dev.off
#' @importFrom lubridate is.POSIXct
#' @importFrom magrittr %>%
#' @importFrom pheatmap pheatmap
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom pvca pvcaBatchAssess
#' @importFrom purrr pmap negate map
#' @importFrom RColorBrewer brewer.pal brewer.pal.info
#' @importFrom rlang :=
#' @importFrom rlang !!
#' @importFrom rlang !!!
#' @importFrom rlang sym syms
#' @importFrom sva ComBat
#' @importFrom tidyr complete nest unnest
#' @importFrom utils combn 
#' @importFrom scales brewer_pal
#' @importFrom stats as.formula complete.cases cor dist hclust sd
#' @importFrom stats ksmooth loess median 
#' @importFrom stats model.matrix prcomp predict reformulate setNames
#' @importFrom tibble remove_rownames rownames_to_column column_to_rownames
#' @importFrom tools file_ext
#' @importFrom viridis viridis_pal
#' @importFrom wesanderson wes_palettes
#' @importFrom WGCNA plotDendroAndColors standardColors
#' 
#' @docType package
#' @name proBatch
#' @aliases proBatch-package
"_PACKAGE"
if(getRversion() >= "2.15.1")  utils::globalVariables(c( "batch_size", 
                                                         "tipping.points", 
                                                         "min_order_value", 
  "data", "batch_total", "fit",  "mean_fit",
  "CV_total","CV_perBatch", "diff_fit", "diff_medians", "sd",
  "median_global", "median_batch", "diff_norm",
  "mean_global", "mean_batch", "diff_means",
  "dateTime", 
  "same_protein", "batch_the_same", 
  "median_run",
  "Var1", "Var2", "label", "weights", "category",
  "Step", "correlation",
  "."))
NULL


