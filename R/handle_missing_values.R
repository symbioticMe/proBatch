handle_missing_values <- function(data_matrix, warning_message,
                                  fill_the_missing = NULL) {
  if (any(is.na(as.vector(data_matrix)))) {
    warning(warning_message)
    if (!is.null(fill_the_missing)) {
      if (!is.numeric(fill_the_missing)) {
        fill_the_missing <- 0
      } else {
        warning(sprintf("filling missing value with %s", fill_the_missing))
        data_matrix[is.na(data_matrix)] <- fill_the_missing
      }
    } else {
      complete_cases <- complete.cases(data_matrix)
      if (ncol(data_matrix) == nrow(data_matrix)) {
        pre_corr_dim <- nrow(data_matrix)
        if (isSymmetric(data_matrix)) {
          warning("removing rows and columns with missing values, as matrix is square")
          if (all(!complete_cases)) {
            warning("removing rows with all values missing")
            all_missing_rows <- apply(data_matrix, 2, function(x) all(is.na(x)))
            bad_rows <- names(all_missing_rows[all_missing_rows])
            good_rows <- setdiff(colnames(data_matrix), bad_rows)
            data_matrix <- (data_matrix[good_rows, good_rows])
          }
          data_matrix <- data_matrix[complete.cases(data_matrix), complete.cases(data_matrix)]
          after_corr_dim <- nrow(data_matrix)
          rem_rows <- pre_corr_dim - after_corr_dim
          warning(sprintf("removed %s rows of the matrix with missing values", rem_rows))
        } else {
          warning("matrix is square, but not symmetric, coincidence or error?")
          data_matrix <- data_matrix[complete.cases(data_matrix), ]
        }
      } else {
        warning("filling value is NULL, removing rows with missing values from data matrix")
        data_matrix <- data_matrix[complete.cases(data_matrix), ]
      }
    }
  }
  return(data_matrix)
}
