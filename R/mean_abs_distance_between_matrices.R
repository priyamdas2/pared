#' Mean Absolute Distance Between Matrices
#'
#' Internal helper function for computing the mean absolute pairwise distance
#' among matrices in a list.
#'
#' @param matrix_list A list of matrices with identical dimensions.
#'
#' @return A numeric value giving the mean pairwise absolute distance.
#'
#' @keywords internal
#' @noRd
mean_abs_distance_between_matrices <- function(matrix_list) {
  
  dims <- sapply(matrix_list, dim)
  
  if (!all(apply(dims, 1, function(x) length(unique(x)) == 1))) {
    stop("All matrices must have the same dimensions.")
  }
  
  K <- length(matrix_list)
  
  if (K < 2) {
    return(0)
  }
  
  pair_indices <- combn(K, 2)
  
  total_dists <- vapply(
    seq_len(ncol(pair_indices)),
    function(k) {
      i <- pair_indices[1, k]
      j <- pair_indices[2, k]
      sum(abs(matrix_list[[i]] - matrix_list[[j]]))
    },
    numeric(1)
  )
  
  mean(total_dists)
}
