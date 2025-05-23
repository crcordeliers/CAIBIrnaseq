#' Get Empty Rows in a Matrix
#'
#' This function identifies the rows in a matrix that have a sum of zero across a specified set of columns.
#'
#' @param matrix A numeric matrix where rows represent observations and columns represent variables.
#' @param colToCheck An integer vector specifying which columns to check for zero sums. Default is 1:8.
#'
#' @return An integer vector of row indices where the sum across the specified columns is zero.
#' @export
getEmptyRows <- function(matrix, colToCheck = 1:8) {
  # Input validation
  if (!is.matrix(matrix) || !is.numeric(matrix)) {
    stop("`matrix` must be a numeric matrix.")
  }

  if (any(colToCheck < 1 | colToCheck > ncol(matrix))) {
    stop("`colToCheck` contains invalid column indices.")
  }

  rowToRm <- which(rowSums(matrix[, colToCheck, drop = FALSE]) == 0)
  return(rowToRm)
}
