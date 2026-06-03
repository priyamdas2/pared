#' GP-Based Pareto Front for Fused LASSO Regression
#'
#' Performs Gaussian-process-based Pareto optimization for Fused LASSO regression
#' and returns a 3D Pareto-front plot. The tuning parameters are the LASSO
#' penalty \eqn{\lambda_1} and the fusion penalty \eqn{\lambda_2}, represented
#' internally as \eqn{\log_{10}(\lambda_1)} and
#' \eqn{\log_{10}(\lambda_2)}. The default objectives are the number of non-zero
#' coefficients, residual sum of squares, and coefficient roughness.
#'
#' @param X Numeric predictor matrix of dimension \eqn{n \times p}.
#' @param y Numeric response vector of length \eqn{n}.
#' @param lb Numeric vector of lower bounds for \eqn{\log_{10}(\lambda_1)}
#'   and \eqn{\log_{10}(\lambda_2)}. Default is \code{c(-3, -3)}.
#' @param ub Numeric vector of upper bounds for \eqn{\log_{10}(\lambda_1)}
#'   and \eqn{\log_{10}(\lambda_2)}. Default is \code{c(1, 1)}.
#' @param Pareto_budget Number of function evaluations used in the Pareto
#'   optimization. Default is \code{100}.
#' @param plot_title Character string giving the plot title. Default is
#'   \code{"Fused LASSO"}.
#' @param plot_marker_color Marker fill color for the Pareto-front points.
#'   Default is \code{'rgba(255, 0, 0, 0.7)'}.
#' @param plot_marker_border_color Marker border color. Default is
#'   \code{'black'}.
#' @param plot_marker_size Marker size. Default is \code{10}.
#' @param plot_marker_border_width Marker border width. Default is \code{2}.
#' @param plot_marker_symbol Marker symbol, such as \code{'circle'},
#'   \code{'square'}, or \code{'diamond'}. Default is \code{'circle'}.
#' @param fontsize_x Font size for the x-axis title. Default is \code{20}.
#' @param fontsize_y Font size for the y-axis title. Default is \code{20}.
#' @param fontsize_z Font size for the z-axis title. Default is \code{20}.
#' @param plot_title_size Font size for the plot title. Default is \code{40}.
#' @param title_pos_x Horizontal plot-title position on the Plotly 0--1 scale.
#'   Default is \code{0.5}.
#' @param title_pos_y Vertical plot-title position on the Plotly 0--1 scale.
#'   Default is \code{0.95}.
#' @param draw_projection Numeric flag, either \code{0} or \code{1}, indicating
#'   whether projection lines to the xy-plane are drawn. Default is \code{1}.
#' @param proj_line_color Projection-line color. Default is \code{'black'}.
#' @param proj_line_width Projection-line width. Default is \code{2}.
#' @param proj_marker_size Marker size for projected points. Default is
#'   \code{2}.
#' @param proj_marker_color Marker color for projected points. Default is
#'   \code{'rgba(255, 0, 0, 0.7)'}.
#' @param proj_marker_border_width Border width for projected points. Default is
#'   \code{2}.
#' @param proj_marker_border_color Border color for projected points. Default is
#'   \code{'black'}.
#' @param objective_directions Character vector specifying whether each objective
#'   should be minimized or maximized. Default is
#'   \code{c("min", "min", "min")}, corresponding to the number of non-zero
#'   coefficients, RSS, and roughness.
#' @param control Optional list of additional arguments passed to
#'   \code{\link{pared_optimize}}. Default is \code{list()}.
#' @param optim_control Optional list of additional arguments passed to
#'   \code{\link[stats]{optim}} when estimating the Fused LASSO coefficients.
#'   Default is \code{list()}.
#' @param round_digits Number of digits used when rounding the returned summary
#'   table and hover labels. Default is \code{3}.
#' @param verbose Logical; if \code{TRUE}, progress and diagnostic messages may
#'   be printed. Default is \code{TRUE}.
#'
#' @return A list with three components:
#' \itemize{
#'   \item \code{figure}: a \code{plotly} 3D scatter plot of the Pareto front;
#'   \item \code{summary_table}: a data frame containing Pareto-optimal tuning
#'   values and objective values;
#'   \item \code{pared_result}: the full object returned by
#'   \code{\link{pared_optimize}}.
#' }
#'
#' @details
#' For each candidate pair \eqn{(\lambda_1,\lambda_2)}, the function estimates
#' the Fused LASSO coefficient vector by minimizing
#' \deqn{
#' \mathcal{L}(\beta) = \frac{1}{2n} \|y - X\beta\|_2^2 + \lambda_1 \|\beta\|_1 + \lambda_2 \sum_{j=2}^{p} |\beta_j - \beta_{j-1}|
#' }
#' The three objectives are the number of non-zero coefficients, residual sum
#' of squares, and roughness \eqn{\sum_j |\beta_{j+1}-\beta_j|}. By default, all
#' three objectives are minimized. The resulting Pareto-optimal solutions are
#' visualized using \code{\link{plot_pared_3d}}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 80
#' p <- 10
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' beta <- c(rep(2, 3), rep(0, 4), rep(-1.5, 3))
#' y <- as.numeric(X %*% beta + rnorm(n))
#'
#' fit <- pared_FLasso(X, y, Pareto_budget = 20)
#' fit$summary_table
#'
#' fit$figure
#' }
#'
#' @seealso \code{\link{pared_optimize}}, \code{\link{plot_pared_3d}}
#'
#' @export


pared_FLasso <- function(X, y,
                         lb = c(-3, -3),
                         ub = c(1, 1),
                         Pareto_budget = 100,
                         plot_title = "Fused LASSO",
                         plot_marker_color = 'rgba(255, 0, 0, 0.7)',
                         plot_marker_border_color = 'black',
                         plot_marker_size = 10,
                         plot_marker_border_width = 2,
                         plot_marker_symbol = 'circle',
                         fontsize_x = 20,
                         fontsize_y = 20,
                         fontsize_z = 20,
                         plot_title_size = 40,
                         title_pos_x = 0.5,
                         title_pos_y = 0.95,
                         draw_projection = 1,
                         proj_line_color = 'black',
                         proj_line_width = 2,
                         proj_marker_size = 2,
                         proj_marker_color = 'rgba(255, 0, 0, 0.7)',
                         proj_marker_border_width = 2,
                         proj_marker_border_color = 'black',
                         objective_directions = c("min", "min", "min"),
                         control = list(),
                         optim_control = list(),
                         round_digits = 3,
                         verbose = TRUE) {

  ###########################################################################
  ## Basic checks
  ###########################################################################

  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }

  y <- as.numeric(y)

  if (nrow(X) != length(y)) {
    stop("nrow(X) must equal length(y).")
  }

  if (length(lb) != 2 || length(ub) != 2) {
    stop("For fused lasso, lb and ub must be length-2 vectors for log10(lambda1), log10(lambda2).")
  }


  ###########################################################################
  ## Local fused-lasso objective in beta
  ###########################################################################

  FLasso_objective_local <- function(X, y, beta, log_lam1, log_lam2) {

    lambda1 <- 10^log_lam1
    lambda2 <- 10^log_lam2

    n <- length(y)

    residuals <- y - X %*% beta
    rss <- sum(residuals^2) / (2 * n)

    l1_penalty <- lambda1 * sum(abs(beta))

    fused_penalty <- lambda2 * sum(abs(diff(beta)))

    rss + l1_penalty + fused_penalty
  }

  ###########################################################################
  ## Objective function for pared_optimize()
  ###########################################################################

  objective_fun <- function(log_lambdas) {

    if (is.null(dim(log_lambdas))) {
      log_lambdas <- matrix(log_lambdas, nrow = 1)
    }

    log_lam1 <- log_lambdas[1]
    log_lam2 <- log_lambdas[2]

    p <- ncol(X)

    beta_objective <- function(beta) {
      FLasso_objective_local(
        X = X,
        y = y,
        beta = beta,
        log_lam1 = log_lam1,
        log_lam2 = log_lam2
      )
    }

    optim_args <- c(
      list(
        par = rep(0, p),
        fn = beta_objective,
        method = "L-BFGS-B",
        lower = rep(-10^3, p),
        upper = rep(10^3, p),
        control = list(maxit = 500)
      ),
      optim_control
    )

    optim_fit <- do.call(stats::optim, optim_args)
    beta_hat <- as.numeric(optim_fit$par)

    num_nonzero <- sum(abs(beta_hat) > 10^(-3))

    residuals <- y - X %*% beta_hat
    rss <- sum(residuals^2) / length(y)

    roughness <- sum(abs(diff(beta_hat)))

    c(num_nonzero, rss, roughness)
  }

  ###########################################################################
  ## Run generic Pareto optimizer
  ###########################################################################

  res <- pared_optimize(
    objective_fun = objective_fun,
    lower = lb,
    upper = ub,
    budget = Pareto_budget,
    parameter_names = c("log10_lambda1", "log10_lambda2"),
    objective_names = c("Num. non-zero coeffs.", "RSS", "Roughness"),
    objective_directions = objective_directions,
    control = control,
    round_digits = round_digits,
    verbose = verbose
  )

  ###########################################################################
  ## Match old-style summary table
  ###########################################################################

  lambda1_vals <- round(10^res$par[, 1], round_digits)
  lambda2_vals <- round(10^res$par[, 2], round_digits)

  Summary.Table <- data.frame(
    lambda1 = lambda1_vals,
    lambda2 = lambda2_vals,
    `Num. non-zero coeffs.` = round(res$value[, 1], round_digits),
    RSS = round(res$value[, 2], round_digits),
    Roughness = round(res$value[, 3], round_digits),
    check.names = FALSE
  )

  res$summary_table <- Summary.Table

  ###########################################################################
  ## Generic 3D plot
  ###########################################################################

  fig <- plot_pared_3d(
    pared_result = res,
    x_objective = "Num. non-zero coeffs.",
    y_objective = "RSS",
    z_objective = "Roughness",
    hover_parameters = c("lambda1", "lambda2"),
    plot_title = plot_title,
    xaxis_title = "No. of non-zero coeffs.",
    yaxis_title = "RSS",
    zaxis_title = "Roughness",
    plot_marker_color = plot_marker_color,
    plot_marker_border_color = plot_marker_border_color,
    plot_marker_size = plot_marker_size,
    plot_marker_border_width = plot_marker_border_width,
    plot_marker_symbol = plot_marker_symbol,
    fontsize_x = fontsize_x,
    fontsize_y = fontsize_y,
    fontsize_z = fontsize_z,
    plot_title_size = plot_title_size,
    title_pos_x = title_pos_x,
    title_pos_y = title_pos_y,
    draw_projection = draw_projection,
    proj_line_color = proj_line_color,
    proj_line_width = proj_line_width,
    proj_marker_size = proj_marker_size,
    proj_marker_color = proj_marker_color,
    proj_marker_border_width = proj_marker_border_width,
    proj_marker_border_color = proj_marker_border_color,
    round_digits = round_digits
  )

  return(list(
    figure = fig,
    summary_table = Summary.Table,
    pared_result = res
  ))
}
