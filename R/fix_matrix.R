#' Stabilize a Precision Matrix
#'
#' Internal helper function for scaling a symmetric matrix to improve positive
#' definiteness while preserving unit diagonal entries.
#'
#' @param A A square numeric matrix.
#' @param denom_factor Numeric scaling factor used in row-wise normalization.
#'
#' @return A symmetric numeric matrix with unit diagonal entries.
#'
#' @keywords internal
#' @noRd
fix_matrix <- function(A, denom_factor) {
  
  p <- nrow(A)
  
  for (cur_row in seq_len(p)) {
    
    cur_sum <- sum(abs(A[cur_row, ])) - 1
    
    if (cur_sum != 0) {
      A[cur_row, ] <- A[cur_row, ] / (denom_factor * cur_sum)
    }
    
    A[cur_row, cur_row] <- 1
  }
  
  A <- (A + t(A)) / 2
  
  return(A)
}