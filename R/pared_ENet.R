#' GP-Based Pareto Front for Elastic-Net Regression
#'
#' Performs Gaussian-process-based Pareto optimization for Elastic-Net regression
#' and returns a 3D Pareto-front plot. The tuning parameters are the Elastic-Net
#' mixing parameter \eqn{\alpha} and the regularization parameter \eqn{\lambda},
#' represented internally as \eqn{\log_{10}(\lambda)}. The default objectives are
#' the number of non-zero coefficients, the L2 norm of the fitted coefficient
#' vector, and the model deviance.
#'
#' @param X Numeric predictor matrix of dimension \eqn{n \times p}.
#' @param y Numeric response vector of length \eqn{n}.
#' @param lb Numeric vector of lower bounds for \eqn{\alpha} and
#'   \eqn{\log_{10}(\lambda)}. Default is \code{c(0, -6)}.
#' @param ub Numeric vector of upper bounds for \eqn{\alpha} and
#'   \eqn{\log_{10}(\lambda)}. Default is \code{c(1, 6)}.
#' @param Pareto_budget Number of function evaluations used in the Pareto
#'   optimization. Default is \code{100}.
#' @param plot_title Character string giving the plot title. Default is
#'   \code{"Elastic-Net"}.
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
#'   coefficients, L2 norm, and deviance.
#' @param control Optional list of additional arguments passed to
#'   \code{\link{pared_optimize}}. Default is \code{list()}.
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
#' For each candidate pair \eqn{(\alpha,\lambda)}, the function fits an
#' Elastic-Net model using \code{\link[glmnet]{glmnet}}. The objective function
#' evaluates three quantities: the number of non-zero coefficients, the L2 norm
#' of the coefficient vector, and the deviance. By default, all three objectives
#' are minimized. The resulting Pareto-optimal solutions are visualized using
#' \code{\link{plot_pared_3d}}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 100
#' p <- 10
#' X <- matrix(rnorm(n * p), nrow = n, ncol = p)
#' beta <- c(2, -1.5, rep(0, p - 2))
#' y <- as.numeric(X %*% beta + rnorm(n))
#'
#' fit <- pared_ENet(X, y, Pareto_budget = 30)
#' fit$summary_table
#'
#' fit$figure
#' }
#'
#' @seealso \code{\link{pared_optimize}}, \code{\link{plot_pared_3d}},
#'   \code{\link[glmnet]{glmnet}}
#'
#' @importFrom glmnet glmnet
#' @export


pared_ENet <- function(X, y,
                       lb = c(0, -6),
                       ub = c(1, 6),
                       Pareto_budget = 100,
                       plot_title = "Elastic-Net",
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
                       round_digits = 3,
                       verbose = TRUE) {

  objective_fun <- function(alpha.loglambda) {

    if (is.null(dim(alpha.loglambda))) {
      alpha.loglambda <- matrix(alpha.loglambda, nrow = 1)
    }

    ElasticNet.fit <- glmnet::glmnet(
      x = X,
      y = y,
      alpha = alpha.loglambda[1],
      lambda = 10^(alpha.loglambda[2])
    )

    g1 <- sum(abs(ElasticNet.fit$beta) > 1e-6)

    g2 <- sqrt(sum(ElasticNet.fit$beta^2))

    Dev <- deviance(ElasticNet.fit)
    g3 <- Dev[1]

    if (length(Dev) > 1) {
      message("Deviance is of length more than 1")
    }

    c(g1, g2, g3)
  }

  res <- pared_optimize(
    objective_fun = objective_fun,
    lower = lb,
    upper = ub,
    budget = Pareto_budget,
    parameter_names = c("alpha", "log10_lambda"),
    objective_names = c("Num. non-zero coeffs.", "L2 norm", "Deviance"),
    objective_directions = objective_directions,
    control = control,
    round_digits = round_digits,
    verbose = verbose
  )

  fig <- plot_pared_3d(
    pared_result = res,
    x_objective = "Num. non-zero coeffs.",
    y_objective = "L2 norm",
    z_objective = "Deviance",
    hover_parameters = c("alpha", "log10_lambda"),
    plot_title = plot_title,
    xaxis_title = "No. of non-zero coeffs.",
    yaxis_title = "L2(\u03B2)",
    zaxis_title = "Deviance",
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
    summary_table = res$summary_table,
    pared_result = res
  ))
}
