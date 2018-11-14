#' Plot boxplots to compare various data normalization steps/approaches WARNING:
#' extremely slow for big dataframes
#'
#' @param list_of_dfs list of data frames of format, specified in `plot_boxplot`
#' @param sample_annotation data matrix with 1) `sample_id_col` (this can be
#'   repeated as row names) 2) biological and 3) technical covariates (batches
#'   etc)
#' @param batch_col column in `sample_annotation` that should be used for
#'   batch comparison
#' @param step normalization step (e.g. `Raw` or `Quantile_normalized` or
#'   `qNorm_ComBat`). Useful if consecutive steps are compared in plots. Note
#'   that in plots these are usually ordered alphabetically, so it's worth
#'   naming with numbers, e.g. `1_raw`, `2_quantile`
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' @seealso \code{\link{plot_boxplot}}
boxplot_all_steps <- function(list_of_dfs, sample_annotation, batch_col,
                              step = NULL){
  if(`*`(dim(list_of_dfs[[1]])[1], dim(list_of_dfs[[1]])[2]) * length(list_of_dfs) > 5*10^5){
    warning('Data matrices are huge, be patient, this might take a while (or crash)')
  }
  add_processing_step <- function(i, list_of_dfs, steps) {
    df = list_of_dfs[[i]]
    df$step = steps[i]
    list_of_dfs[[i]] = df
  }
  if (!is.null(step) |
      (any(sapply(list_of_dfs, function(df) (!'step' %in% names(df)))))){
    list_of_dfs = lapply(1:length(list_of_dfs), add_processing_step, list_of_dfs, step)
  }
  
  joined_proteome = do.call(rbind, list_of_dfs)
  gg = gg_boxplot(joined_proteome, sample_annotation, batch_col) +
    facet_grid(.~step)
  return(gg)
}