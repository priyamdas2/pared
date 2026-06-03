#' General GP-Based Pareto Optimization
#'
#' Performs Gaussian-process-based Pareto optimization for a user-defined
#' multi-objective function. The user supplies an objective function, lower and
#' upper bounds for the search parameters, optional parameter and objective
#' names, and optimization directions for each objective.
#'
#' @param objective_fun A user-defined function taking a numeric parameter
#'   vector as input and returning a numeric vector of objective values.
#' @param lower Numeric vector of lower bounds for the search parameters.
#' @param upper Numeric vector of upper bounds for the search parameters.
#' @param budget Number of function evaluations used in the Pareto optimization.
#'   Default is \code{100}.
#' @param parameter_names Character vector giving names for the search
#'   parameters. Default is \code{NULL}; if \code{NULL}, names are assigned as
#'   \code{paste0("par", seq_len(length(lower)))}.
#' @param objective_names Character vector giving names for the objectives.
#'   Default is \code{NULL}; if \code{NULL}, names are assigned as
#'   \code{paste0("objective", seq_len(q))}, where \code{q} is the number of
#'   objective values returned by \code{objective_fun}.
#' @param objective_directions Character vector specifying whether each
#'   objective should be minimized or maximized. Each entry must be either
#'   \code{"min"} or \code{"max"}. Default is \code{NULL}; if \code{NULL},
#'   all objectives are minimized, i.e., \code{rep("min", q)}.
#' @param control Optional list of additional arguments passed to
#'   \code{\link[GPareto]{easyGParetoptim}}. Default is \code{list()}.
#' @param round_digits Number of digits used when rounding the returned summary
#'   table and diagnostic messages. Default is \code{3}.
#' @param remove_duplicate_objectives Logical; if \code{TRUE}, duplicate
#'   objective rows are removed from the returned Pareto summary. Default is
#'   \code{TRUE}.
#' @param return_raw Logical; if \code{TRUE}, raw parameter and objective
#'   matrices from the optimizer are included in the output. Default is
#'   \code{TRUE}.
#' @param verbose Logical; if \code{TRUE}, diagnostic messages are printed when
#'   the objective function fails or returns non-finite values. Default is
#'   \code{TRUE}.
#' @param ... Additional arguments passed to \code{objective_fun}.
#'
#' @return An object of class \code{"pared_result"}, returned as a list with
#' the following components:
#' \itemize{
#'   \item \code{par}: matrix of Pareto-optimal parameter values after optional
#'   duplicate-objective removal;
#'   \item \code{value}: matrix of objective values on the natural display
#'   scale;
#'   \item \code{value_for_optimization}: matrix of objective values used
#'   internally for optimization, with maximization objectives sign-reversed;
#'   \item \code{summary_table}: data frame combining rounded parameter values
#'   and rounded objective values;
#'   \item \code{parameter_names}: parameter names used in the output;
#'   \item \code{objective_names}: objective names used in the output;
#'   \item \code{objective_directions}: optimization directions used for the
#'   objectives;
#'   \item \code{lower}: lower search bounds used in the optimization;
#'   \item \code{upper}: upper search bounds used in the optimization;
#'   \item \code{budget}: optimization budget used in the run;
#'   \item \code{control}: control list passed to the optimizer;
#'   \item \code{gp_result}: raw object returned by
#'   \code{\link[GPareto]{easyGParetoptim}}.
#' }
#' If \code{return_raw = TRUE}, the output also includes:
#' \itemize{
#'   \item \code{par_raw}: raw parameter matrix returned by the optimizer;
#'   \item \code{value_for_optimization_raw}: raw objective matrix used
#'   internally for optimization;
#'   \item \code{value_raw}: raw objective matrix converted back to the natural
#'   display scale.
#' }
#'
#' @details
#' The function first evaluates \code{objective_fun} at the midpoint of the
#' search space to determine the number of objectives and to check that the
#' objective function returns finite numeric values. Objectives marked as
#' \code{"max"} are internally multiplied by \code{-1}, so that all objectives
#' are passed to \code{\link[GPareto]{easyGParetoptim}} as minimization
#' objectives. The returned \code{value}, \code{value_raw}, and
#' \code{summary_table} are converted back to the original objective scale.
#'
#' If \code{objective_fun} fails or returns non-finite values at a candidate
#' point during optimization, large penalty values are returned for that point.
#' This allows the optimization to continue while discouraging infeasible or
#' numerically unstable regions of the search space.
#'
#' @examples
#' \dontrun{
#' toy_objective <- function(z) c(Mean = mean(z),
#'                               Harmonic_mean = length(z) / sum(1 / z),
#'                               MAD_from_1.5 = mean(abs(z - 1.5)))
#'
#' set.seed(1)
#' toy_res <- pared_optimize(objective_fun = toy_objective,
#'                           lower = c(1, 1, 1),
#'                           upper = c(2, 2, 2),
#'                           budget = 30,
#'                           parameter_names = c("z1", "z2", "z3"),
#'                           objective_names = c("Mean", "Harmonic mean",
#'                                               "MAD around 1.5"),
#'                           objective_directions = c("max", "min", "min"),
#'                           verbose = FALSE)
#'
#' print(toy_res$summary_table)
#' 
#' 
#' plot_pared_2d(
#'  toy_res,
#'  x_objective = "Mean",
#'  y_objective = "MAD around 1.5",
#'  plot_title = "Toy Pareto front"
#')
#' 
#' plot_pared_3d(toy_res,
#'               x_objective = "Mean",
#'               y_objective = "Harmonic mean",
#'               z_objective = "MAD around 1.5",
#'               plot_title = "Toy Pareto front",
#'               plot_marker_color = "blue"
#'               )
#' }
#'
#' @seealso \code{\link{plot_pared_2d}}, \code{\link{plot_pared_3d}},
#'   \code{\link[GPareto]{easyGParetoptim}}
#'
#' @export


pared_optimize <- function(objective_fun,
                           lower,
                           upper,
                           budget = 100,
                           parameter_names = NULL,
                           objective_names = NULL,
                           objective_directions = NULL,
                           control = list(),
                           round_digits = 3,
                           remove_duplicate_objectives = TRUE,
                           return_raw = TRUE,
                           verbose = TRUE,
                           ...) {
  
  ###########################################################################
  ## Basic checks
  ###########################################################################
  
  if (!is.function(objective_fun)) {
    stop("objective_fun must be a function.")
  }
  
  if (!is.numeric(lower) || !is.numeric(upper)) {
    stop("lower and upper must be numeric vectors.")
  }
  
  if (length(lower) != length(upper)) {
    stop("lower and upper must have the same length.")
  }
  
  if (any(!is.finite(lower)) || any(!is.finite(upper))) {
    stop("lower and upper must contain only finite values.")
  }
  
  if (any(lower >= upper)) {
    stop("Each element of lower must be strictly smaller than the corresponding element of upper.")
  }
  
  if (!is.numeric(budget) || length(budget) != 1 || budget <= 0) {
    stop("budget must be a positive numeric scalar.")
  }
  
  budget <- as.integer(budget)
  
  d <- length(lower)
  
  if (is.null(parameter_names)) {
    parameter_names <- paste0("par", seq_len(d))
  }
  
  if (length(parameter_names) != d) {
    stop("parameter_names must have the same length as lower and upper.")
  }
  
  if (!is.list(control)) {
    stop("control must be a list.")
  }
  
  ###########################################################################
  ## First evaluation to infer number of objectives
  ###########################################################################
  
  test_par <- (lower + upper) / 2
  
  test_val <- tryCatch(
    {
      objective_fun(test_par, ...)
    },
    error = function(e) {
      stop(
        paste0(
          "objective_fun failed at the midpoint of the search space. ",
          "Original error: ", e$message
        )
      )
    }
  )
  
  test_val <- as.numeric(test_val)
  
  if (length(test_val) == 0) {
    stop("objective_fun must return at least one numeric objective value.")
  }
  
  if (any(!is.finite(test_val))) {
    stop("objective_fun must return finite numeric objective values at the midpoint of the search space.")
  }
  
  q <- length(test_val)
  
  if (is.null(objective_names)) {
    objective_names <- paste0("objective", seq_len(q))
  }
  
  if (length(objective_names) != q) {
    stop("objective_names must have the same length as the output of objective_fun.")
  }
  
  if (is.null(objective_directions)) {
    objective_directions <- rep("min", q)
  }
  
  objective_directions <- tolower(objective_directions)
  
  if (length(objective_directions) != q) {
    stop("objective_directions must have the same length as the output of objective_fun.")
  }
  
  if (!all(objective_directions %in% c("min", "max"))) {
    stop("objective_directions must contain only 'min' or 'max'.")
  }
  
  maximize_index <- which(objective_directions == "max")
  
  ###########################################################################
  ## Internal wrapper for GPareto
  ###########################################################################
  
  wrapped_objective_fun <- function(par) {
    
    penalty_vals_for_optimization <- rep(.Machine$double.xmax ^ 0.25, q)
    
    vals <- tryCatch(
      {
        objective_fun(par, ...)
      },
      error = function(e) {
        if (verbose) {
          message("objective_fun failed at par = ",
                  paste(round(par, round_digits), collapse = ", "),
                  ". Returning large penalty values.")
        }
        return(NULL)
      }
    )
    
    if (is.null(vals)) {
      penalty_mat <- matrix(penalty_vals_for_optimization, nrow = 1)
      dimnames(penalty_mat) <- NULL
      return(penalty_mat)
    }
    
    vals <- as.numeric(vals)
    
    if (length(vals) != q) {
      stop(
        paste0(
          "objective_fun returned ", length(vals),
          " values, but it must always return ", q, " values."
        )
      )
    }
    
    if (any(!is.finite(vals))) {
      if (verbose) {
        message("objective_fun returned non-finite values at par = ",
                paste(round(par, round_digits), collapse = ", "),
                ". Returning large penalty values.")
      }
      penalty_mat <- matrix(penalty_vals_for_optimization, nrow = 1)
      dimnames(penalty_mat) <- NULL
      return(penalty_mat)
    }
    
    vals_for_optimization <- vals
    
    if (length(maximize_index) > 0) {
      vals_for_optimization[maximize_index] <- -vals_for_optimization[maximize_index]
    }
    
    vals_for_optimization <- matrix(vals_for_optimization, nrow = 1)
    dimnames(vals_for_optimization) <- NULL
    
    return(vals_for_optimization)
  }
  
  ###########################################################################
  ## Run GP-based Pareto optimization
  ###########################################################################
  
  gp_args <- c(
    list(
      fn = wrapped_objective_fun,
      budget = budget,
      lower = lower,
      upper = upper
    ),
    control
  )
  
  gp_fit <- do.call(GPareto::easyGParetoptim, gp_args)
  
  ###########################################################################
  ## Extract parameters and objective values
  ###########################################################################
  
  par_raw <- gp_fit$par
  value_for_optimization_raw <- gp_fit$value
  
  if (is.null(dim(par_raw))) {
    par_raw <- matrix(par_raw, nrow = 1)
  }
  
  if (is.null(dim(value_for_optimization_raw))) {
    value_for_optimization_raw <- matrix(value_for_optimization_raw, nrow = 1)
  }
  
  ###########################################################################
  ## Remove duplicate objective rows, if requested
  ###########################################################################
  
  if (remove_duplicate_objectives) {
    unique_rows <- !duplicated(value_for_optimization_raw)
  } else {
    unique_rows <- rep(TRUE, nrow(value_for_optimization_raw))
  }
  
  par_unique <- par_raw[unique_rows, , drop = FALSE]
  value_for_optimization_unique <- value_for_optimization_raw[unique_rows, , drop = FALSE]
  
  ###########################################################################
  ## Convert maximization objectives back to natural display scale
  ###########################################################################
  
  value_display <- value_for_optimization_unique
  
  if (length(maximize_index) > 0) {
    value_display[, maximize_index] <- -value_display[, maximize_index]
  }
  
  ###########################################################################
  ## Construct summary table
  ###########################################################################
  
  colnames(par_unique) <- parameter_names
  colnames(value_display) <- objective_names
  colnames(value_for_optimization_unique) <- objective_names
  
  summary_table <- data.frame(
    round(par_unique, round_digits),
    round(value_display, round_digits),
    check.names = FALSE
  )
  
  ###########################################################################
  ## Return object
  ###########################################################################
  
  out <- list(
    par = par_unique,
    value = value_display,
    value_for_optimization = value_for_optimization_unique,
    summary_table = summary_table,
    parameter_names = parameter_names,
    objective_names = objective_names,
    objective_directions = objective_directions,
    lower = lower,
    upper = upper,
    budget = budget,
    control = control,
    gp_result = gp_fit
  )
  
  if (return_raw) {
    out$par_raw <- par_raw
    out$value_for_optimization_raw <- value_for_optimization_raw
    
    value_raw_display <- value_for_optimization_raw
    if (length(maximize_index) > 0) {
      value_raw_display[, maximize_index] <- -value_raw_display[, maximize_index]
    }
    colnames(value_raw_display) <- objective_names
    out$value_raw <- value_raw_display
  }
  
  class(out) <- "pared_result"
  
  return(out)
}