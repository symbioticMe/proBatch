#'proBatch: A package for diagnostics and correction of batch effects, primarily
#'in proteomics
#'
#'It adresses the following needs: \itemize{ \item prepare the original data
#'(e.g. OpenSWATH output matrix and sample annotation file) for analysis.
#'However, you might need to use `SWATH2stats` additionally \item Diagnose batch
#'effects, sample-wide and feature-level \item Correct for batch effects
#'(normalize the data). Other useful package for this purpose is `Normalyzer`. }
#'
#'#' To learn more about proBatch, start with the vignettes:
#'`browseVignettes(package = "proBatch")`
#'
#'@section common arguments to the functions:
#'@param df_long data frame where each row is a single feature in a single
#'  sample, thus it has minimally, `sample_id_col`, `feature_id_column` and
#'  `measure_column`, but usually also `m_score` (in OpenSWATH output result
#'  file)
#'@param data_matrix features (in rows) vs samples (in columns) matrix, with
#'  feature IDs in rownames and file/sample names as colnames. Usually the log
#'  transformed version of the original data
#'@param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'  repeated as row names) 2) biological and 3) technical covariates (batches
#'  etc)
#'@param sample_id_col name of the column in sample_annotation file, where the
#'  filenames (colnames of the data matrix are found)
#'@param batch_column column in `sample_annotation` that should be used for
#'  batch comparison
#'@param measure_column if `df_long` is among the parameters, it is the column
#'  with expression/abundance/intensity; otherwise, it is used internally for
#'  consistency
#'
#'@docType package
#'@name proBatch
NULL
