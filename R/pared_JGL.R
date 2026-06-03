#' GP-Based Pareto Front for Joint Graphical LASSO
#'
#' Performs Gaussian-process-based Pareto optimization for Joint Graphical LASSO
#' tuning and returns a 3D Pareto-front plot. The tuning parameters are
#' \eqn{\lambda_1} and \eqn{\lambda_2}, represented internally as
#' \eqn{\log_{10}(\lambda_1)} and \eqn{\log_{10}(\lambda_2)}. The function
#' supports both group graphical LASSO and fused graphical LASSO penalties.
#'
#' @param sample_list A list of numeric data matrices, one for each class or
#'   group. Each matrix should have observations in rows and variables in
#'   columns. All matrices must have the same number of columns.
#' @param lb Numeric vector of lower bounds for \eqn{\log_{10}(\lambda_1)}
#'   and \eqn{\log_{10}(\lambda_2)}. Default is \code{c(-2, -2)}.
#' @param ub Numeric vector of upper bounds for \eqn{\log_{10}(\lambda_1)}
#'   and \eqn{\log_{10}(\lambda_2)}. Default is \code{c(0, 0)}.
#' @param Pareto_budget Number of function evaluations used in the Pareto
#'   optimization. Default is \code{100}.
#' @param method Character string specifying the Joint Graphical LASSO penalty.
#'   Must be either \code{"group"} or \code{"fused"}. Default is
#'   \code{"group"}.
#' @param plot_title Character string giving the plot title. Default is
#'   \code{NULL}; if \code{NULL}, the title is set to
#'   \code{"JGL: Group LASSO"} when \code{method = "group"} and
#'   \code{"JGL: Fused LASSO"} when \code{method = "fused"}.
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
#' @param control Optional list of additional arguments passed to
#'   \code{\link{pared_optimize}}. Default is \code{list()}.
#' @param jgl_control Optional list of additional arguments passed to the
#'   underlying \code{JGL} fitting routine. Default is \code{list()}.
#' @param round_digits Number of digits used when rounding the returned summary
#'   table and hover labels. Default is \code{3}.
#' @param verbose Logical; if \code{TRUE}, progress and diagnostic messages may
#'   be printed. Default is \code{FALSE}.
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
#' For each candidate pair \eqn{(\lambda_1,\lambda_2)}, the function fits a
#' Joint Graphical LASSO model using the selected penalty type. When
#' \code{method = "group"}, the three objectives are total number of edges,
#' number of shared edges, and AIC; the corresponding optimization directions
#' are minimize, maximize, and minimize. When \code{method = "fused"}, the three
#' objectives are total number of edges, mean absolute deviation among the
#' estimated precision matrices, and AIC; all three are minimized. The resulting
#' Pareto-optimal solutions are visualized using \code{\link{plot_pared_3d}}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 60
#' p <- 5
#' sample_list <- list(
#'   matrix(rnorm(n * p), nrow = n, ncol = p),
#'   matrix(rnorm(n * p), nrow = n, ncol = p)
#' )
#'
#' fit <- pared_JGL(sample_list, Pareto_budget = 20, method = "group")
#' fit$summary_table
#'
#' fit$figure
#' }
#'
#' @seealso \code{\link{pared_optimize}}, \code{\link{plot_pared_3d}}
#'
#' @export


pared_JGL <- function(sample_list,
                      lb = c(-2, -2),
                      ub = c(0, 0),
                      Pareto_budget = 100,
                      method = "group",
                      plot_title = NULL,
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
                      control = list(),
                      jgl_control = list(),
                      round_digits = 3,
                      verbose = FALSE) {
  
  if (!method %in% c("group", "fused")) {
    stop("method must be either 'group' or 'fused'.")
  }
  
  if (is.null(plot_title)) {
    if (method == "group") {
      plot_title <- "JGL: Group LASSO"
    } else {
      plot_title <- "JGL: Fused LASSO"
    }
  }
  
  C <- length(sample_list)
  p <- ncol(sample_list[[1]])
  
  mean_abs_distance_between_matrices_local <- function(matrix_list) {
    
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
  
  compute_jgl_aic <- function(Precision.est, sample_list) {
    
    C <- length(sample_list)
    p <- ncol(sample_list[[1]])
    
    AIC.sum <- 0
    num.non.zero.edges <- rep(0, C)
    non.zero.edge.serials <- vector("list", C)
    
    for (c in seq_len(C)) {
      
      sample.size <- nrow(sample_list[[c]])
      
      temp <- which(abs(Precision.est[[c]]) > 0.0001)
      non.zero.edge.serials[[c]] <- temp
      
      num.non.zero.edges[c] <- (length(temp) - p) / 2
      
      S.mat <- t(sample_list[[c]]) %*% sample_list[[c]]
      
      LL <- - (sample.size / 2) * (
        log(det(Precision.est[[c]])) -
          psych::tr(S.mat %*% Precision.est[[c]])
      )
      
      AIC <- 2 * num.non.zero.edges[c] - 2 * LL
      
      AIC.sum <- AIC.sum + AIC
    }
    
    common.elements <- non.zero.edge.serials[[1]]
    
    if (C > 1) {
      for (c in 2:C) {
        common.elements <- intersect(common.elements, non.zero.edge.serials[[c]])
      }
    }
    
    num.common.edges <- (length(common.elements) - p) / 2
    
    list(
      total_edges = sum(num.non.zero.edges),
      shared_edges = num.common.edges,
      AIC = AIC.sum
    )
  }
  
  objective_fun <- function(log_lambdas) {
    
    if (is.null(dim(log_lambdas))) {
      log_lambdas <- matrix(log_lambdas, nrow = 1)
    }
    
    lambda1 <- 10^log_lambdas[1]
    lambda2 <- 10^log_lambdas[2]
    
    jgl_args <- c(
      list(
        Y = sample_list,
        penalty = method,
        lambda1 = lambda1,
        lambda2 = lambda2
      ),
      jgl_control
    )
    
    JGL.fit <- do.call(JGL, jgl_args)
    
    Precision.est <- JGL.fit$theta
    
    num.nodes.now <- sqrt(length(Precision.est[[1]]))
    
    if (num.nodes.now < p) {
      if (verbose) {
        message(sprintf(
          "At (lambda1, lambda2) = (%.3f, %.3f), the number of estimated nodes is less than p. Returning penalty values.",
          lambda1, lambda2
        ))
      }
      
      if (method == "group") {
        return(c(1e10, -1e10, 1e10))
      } else {
        return(c(1e10, 1e10, 1e10))
      }
    }
    
    aic_info <- compute_jgl_aic(
      Precision.est = Precision.est,
      sample_list = sample_list
    )
    
    total_edges <- aic_info$total_edges
    shared_edges <- aic_info$shared_edges
    AIC <- aic_info$AIC
    
    if (method == "group") {
      out <- c(total_edges, shared_edges, AIC)
    } else {
      precision_MAD <- mean_abs_distance_between_matrices_local(Precision.est)
      out <- c(total_edges, precision_MAD, AIC)
    }
    
    return(out)
  }
  
  if (method == "group") {
    objective_names <- c("total edges", "shared edges", "AIC")
    objective_directions <- c("min", "max", "min")
  } else {
    objective_names <- c("total edges", "precision MAD", "AIC")
    objective_directions <- c("min", "min", "min")
  }
  
  res <- pared_optimize(
    objective_fun = objective_fun,
    lower = lb,
    upper = ub,
    budget = Pareto_budget,
    parameter_names = c("log10_lambda1", "log10_lambda2"),
    objective_names = objective_names,
    objective_directions = objective_directions,
    control = control,
    round_digits = round_digits,
    verbose = verbose
  )
  
  lambda1_vals <- round(10^res$par[, 1], round_digits)
  lambda2_vals <- round(10^res$par[, 2], round_digits)
  
  if (method == "group") {
    
    Summary.Table <- data.frame(
      lambda1 = lambda1_vals,
      lambda2 = lambda2_vals,
      `total edges` = round(res$value[, 1], round_digits),
      `shared edges` = round(res$value[, 2], round_digits),
      AIC = round(res$value[, 3], round_digits),
      check.names = FALSE
    )
    
    res$summary_table <- Summary.Table
    
    fig <- plot_pared_3d(
      pared_result = res,
      x_objective = "total edges",
      y_objective = "shared edges",
      z_objective = "AIC",
      hover_parameters = c("lambda1", "lambda2"),
      plot_title = plot_title,
      xaxis_title = "Total no. of edges",
      yaxis_title = "Shared no. of edges",
      zaxis_title = "AIC",
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
    
  } else {
    
    Summary.Table <- data.frame(
      lambda1 = lambda1_vals,
      lambda2 = lambda2_vals,
      `total edges` = round(res$value[, 1], round_digits),
      `precision MAD` = round(res$value[, 2], round_digits),
      AIC = round(res$value[, 3], round_digits),
      check.names = FALSE
    )
    
    res$summary_table <- Summary.Table
    
    fig <- plot_pared_3d(
      pared_result = res,
      x_objective = "total edges",
      y_objective = "precision MAD",
      z_objective = "AIC",
      hover_parameters = c("lambda1", "lambda2"),
      plot_title = plot_title,
      xaxis_title = "Total no. of edges",
      yaxis_title = "Precision MAD",
      zaxis_title = "AIC",
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
  }
  
  return(list(
    figure = fig,
    summary_table = Summary.Table,
    pared_result = res
  ))
}