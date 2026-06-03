#' Generate Synthetic Gaussian Samples with Evolving Graph Structures
#'
#' Generates four synthetic multivariate Gaussian data matrices with evolving
#' precision-matrix graph structures. The dimension is fixed at 20 variables.
#' Starting from a banded AR(2)-type precision matrix, edges are progressively
#' added and removed to create four related group-specific graph structures.
#'
#' @param sample_sizes Numeric vector of length 4 specifying the number of
#'   observations to generate for the four groups. Default is
#'   \code{c(20, 40, 60, 80)}.
#' @param rand_seed Integer random seed for reproducibility. Default is
#'   \code{1}.
#'
#' @return A list of length 4. Each element is a numeric matrix with observations
#'   in rows and variables in columns.
#'
#' @details
#' The function constructs four precision matrices. The first matrix is a
#' banded Toeplitz matrix resembling an AR(2) structure. The remaining matrices
#' are obtained by adding and removing edges from the previous structure. The
#' matrices are adjusted using \code{fix_matrix()} to ensure positive
#' definiteness, and samples are generated from multivariate Gaussian
#' distributions using \code{\link[mvtnorm]{rmvnorm}}.
#'
#' @examples
#' \dontrun{
#' sample_list <- generate_sample(sample_sizes = c(20, 40, 60, 80),
#'                                rand_seed = 1)
#' length(sample_list)
#' dim(sample_list[[1]])
#' }
#'
#' @seealso \code{\link{pared_JGL}}
#'
#' @export

generate_sample <- function(sample_sizes = c(20, 40, 60, 80),
                            rand_seed = 1) {
  
  if (length(sample_sizes) != 4) {
    stop("sample_sizes must be a numeric vector of length 4.")
  }
  
  if (any(sample_sizes <= 0)) {
    stop("All sample sizes must be positive.")
  }
  
  if (!exists("fix_matrix", mode = "function")) {
    stop("fix_matrix() was not found. Please ensure it is available in the package.")
  }
  
  set.seed(rand_seed)
  
  C <- 4
  p <- 20
  
  ## A1 is an AR(2)-type graph.
  A1 <- toeplitz(c(1, 0.5, 0.4, rep(0, p - 3)))
  
  n_edges <- (sum(A1 != 0) - p) / 2
  n_possible <- p * (p - 1) / 2
  
  pos_A1 <- which(upper.tri(A1) & A1 != 0, arr.ind = TRUE)
  zero_A1 <- which(upper.tri(A1) & A1 == 0, arr.ind = TRUE)
  
  pos_inds <- sample(seq_len(nrow(pos_A1)), n_edges)
  zero_inds <- sample(seq_len(nrow(zero_A1)), n_possible - n_edges)
  
  A2 <- A1
  
  ## Add 5 new edges and remove 5 existing edges.
  for (j in 1:5) {
    
    sign_j <- ifelse(runif(1) > 0.5, 1, -1)
    nonzero_val <- runif(1, 0.4, 0.6) * sign_j
    
    r_add <- zero_A1[zero_inds[j], 1]
    c_add <- zero_A1[zero_inds[j], 2]
    A2[r_add, c_add] <- nonzero_val
    A2[c_add, r_add] <- nonzero_val
    
    r_rem <- pos_A1[pos_inds[j], 1]
    c_rem <- pos_A1[pos_inds[j], 2]
    A2[r_rem, c_rem] <- 0
    A2[c_rem, r_rem] <- 0
  }
  
  A3 <- A2
  
  ## Add 10 and remove 10 additional edges.
  for (j in 6:15) {
    
    sign_j <- ifelse(runif(1) > 0.5, 1, -1)
    nonzero_val <- runif(1, 0.4, 0.6) * sign_j
    
    r_add <- zero_A1[zero_inds[j], 1]
    c_add <- zero_A1[zero_inds[j], 2]
    A3[r_add, c_add] <- nonzero_val
    A3[c_add, r_add] <- nonzero_val
    
    r_rem <- pos_A1[pos_inds[j], 1]
    c_rem <- pos_A1[pos_inds[j], 2]
    A3[r_rem, c_rem] <- 0
    A3[c_rem, r_rem] <- 0
  }
  
  A4 <- A3
  
  ## Add 5 and remove 5 additional edges.
  for (j in 16:20) {
    
    sign_j <- ifelse(runif(1) > 0.5, 1, -1)
    nonzero_val <- runif(1, 0.4, 0.6) * sign_j
    
    r_add <- zero_A1[zero_inds[j], 1]
    c_add <- zero_A1[zero_inds[j], 2]
    A4[r_add, c_add] <- nonzero_val
    A4[c_add, r_add] <- nonzero_val
    
    r_rem <- pos_A1[pos_inds[j], 1]
    c_rem <- pos_A1[pos_inds[j], 2]
    A4[r_rem, c_rem] <- 0
    A4[c_rem, r_rem] <- 0
  }
  
  ## Ensure positive definiteness.
  A2 <- fix_matrix(A2, 1)
  A3 <- fix_matrix(A3, 1)
  A4 <- fix_matrix(A4, 1)
  
  Precisions_true <- array(0, dim = c(p, p, C))
  Precisions_true[, , 1] <- A1
  Precisions_true[, , 2] <- A2
  Precisions_true[, , 3] <- A3
  Precisions_true[, , 4] <- A4
  
  sample_list <- vector("list", C)
  
  for (c in seq_len(C)) {
    sample_list[[c]] <- mvtnorm::rmvnorm(
      n = sample_sizes[c],
      mean = rep(0, p),
      sigma = solve(Precisions_true[, , c])
    )
  }
  
  return(sample_list)
}

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