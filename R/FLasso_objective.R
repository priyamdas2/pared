#' Fused LASSO objective helper
#'
#' Exported for advanced users and internal use by pared_FLasso().
#'
#' @keywords internal
#' @export
#' @noRd
FLasso_objective <- function(X, y, beta, log_lam1, log_lam2) {
  
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  
  y <- as.numeric(y)
  beta <- as.numeric(beta)
  
  if (nrow(X) != length(y)) {
    stop("nrow(X) must equal length(y).")
  }
  
  if (ncol(X) != length(beta)) {
    stop("ncol(X) must equal length(beta).")
  }
  
  lambda1 <- 10^log_lam1
  lambda2 <- 10^log_lam2
  
  n <- length(y)
  residuals <- y - X %*% beta
  
  rss <- sum(residuals^2) / (2 * n)
  l1_penalty <- lambda1 * sum(abs(beta))
  fused_penalty <- lambda2 * sum(abs(diff(beta)))
  
  rss + l1_penalty + fused_penalty
}
