handle_missing_values <- function(data_matrix, warning_message, fill_the_missing = NULL) {
  if (any(is.na(as.vector(data_matrix)))){
    warning(warning_message)
    if(!is.null(fill_the_missing)){
      if (!is.numeric(fill_the_missing)){
        fill_the_missing = 0
      } else{
        warning(sprintf('filling missing value with %s', fill_the_missing))
        data_matrix[is.na(data_matrix)] = fill_the_missing
      }
    } else {
      warning('filling value is NULL, removing features with missing values')
      data_matrix = data_matrix[complete.cases(data_matrix),]
    }
  }
  return(data_matrix)
}