#' Two-Dimensional Pareto Front Plot
#'
#' Creates a two-dimensional Plotly scatter plot from a \code{pared_result}
#' object returned by \code{\link{pared_optimize}}. The user selects two
#' objective columns from the Pareto summary table to display on the x- and
#' y-axes.
#'
#' @param pared_result An object of class \code{"pared_result"}, typically
#'   returned by \code{\link{pared_optimize}}.
#' @param x_objective Character string giving the name of the objective column
#'   to display on the x-axis.
#' @param y_objective Character string giving the name of the objective column
#'   to display on the y-axis.
#' @param hover_parameters Character vector giving the parameter columns to
#'   display in the hover text. Default is \code{NULL}, in which case
#'   \code{pared_result$parameter_names} is used.
#' @param plot_title Character string giving the plot title. Default is
#'   \code{"Pareto front"}.
#' @param xaxis_title Character string giving the x-axis title. Default is
#'   \code{NULL}, in which case \code{x_objective} is used.
#' @param yaxis_title Character string giving the y-axis title. Default is
#'   \code{NULL}, in which case \code{y_objective} is used.
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
#' @param axis_line_color Axis-line color. Default is \code{"black"}.
#' @param axis_line_width Axis-line width. Default is \code{2}.
#' @param show_axis_line Logical; if \code{TRUE}, axis lines are shown.
#'   Default is \code{TRUE}.
#' @param mirror_axis_line Logical; if \code{TRUE}, axis lines are mirrored on
#'   the opposite sides of the plotting region. Default is \code{TRUE}.
#' @param show_grid Logical; if \code{TRUE}, grid lines are shown. Default is
#'   \code{TRUE}.
#' @param grid_color Grid-line color. Default is
#'   \code{"rgba(220, 220, 220, 0.7)"}.
#' @param plot_title_size Font size for the plot title. Default is \code{20}.
#' @param title_pos_x Horizontal plot-title position on the Plotly 0--1 scale.
#'   Default is \code{0.5}.
#' @param title_pos_y Vertical plot-title position on the Plotly 0--1 scale.
#'   Default is \code{0.90}.
#' @param round_digits Number of digits used when rounding hover-label values.
#'   Default is \code{3}.
#'
#' @return A \code{plotly} object containing a two-dimensional Pareto-front
#'   scatter plot.
#'
#' @details
#' The function extracts the selected objective columns from
#' \code{pared_result$summary_table}. Hover text is constructed from the
#' selected parameter columns. This function is useful for visualizing pairwise
#' projections of Pareto-optimal solutions, especially when more than two
#' objectives are available.
#'
#' @examples
#' \dontrun{
#' toy_objective <- function(z) c(
#'   Mean = mean(z),
#'   Harmonic_mean = length(z) / sum(1 / z),
#'   MAD_from_1.5 = mean(abs(z - 1.5))
#' )
#'
#' set.seed(1)
#' toy_res <- pared_optimize(
#'   objective_fun = toy_objective,
#'   lower = c(1, 1, 1),
#'   upper = c(2, 2, 2),
#'   budget = 30,
#'   parameter_names = c("z1", "z2", "z3"),
#'   objective_names = c("Mean", "Harmonic mean", "MAD around 1.5"),
#'   objective_directions = c("max", "min", "min"),
#'   verbose = FALSE
#' )
#'
#'
#' toy_res$summary_table
#'
#' plot_pared_2d(
#'  toy_res,
#'  x_objective = "Mean",
#'  y_objective = "MAD around 1.5",
#'  plot_title = "Toy Pareto front"
#')
#'}
#' @seealso \code{\link{pared_optimize}}, \code{\link{plot_pared_3d}}
#'
#' @export


plot_pared_2d <- function(pared_result,
                          x_objective,
                          y_objective,
                          hover_parameters = NULL,
                          plot_title = "Pareto front",
                          xaxis_title = NULL,
                          yaxis_title = NULL,
                          plot_marker_color = 'rgba(255, 0, 0, 0.7)',
                          plot_marker_border_color = 'black',
                          plot_marker_size = 10,
                          plot_marker_border_width = 2,
                          plot_marker_symbol = 'circle',
                          fontsize_x = 20,
                          fontsize_y = 20,
                          axis_line_color = "black",
                          axis_line_width = 2,
                          show_axis_line = TRUE,
                          mirror_axis_line = TRUE,
                          show_grid = TRUE,
                          grid_color = "rgba(220, 220, 220, 0.7)",
                          plot_title_size = 20,
                          title_pos_x = 0.5,
                          title_pos_y = 0.90,
                          round_digits = 3) {

  if (!inherits(pared_result, "pared_result")) {
    stop("pared_result must be an object returned by pared_optimize().")
  }

  summary_table <- pared_result$summary_table

  required_cols <- c(x_objective, y_objective)

  missing_cols <- setdiff(required_cols, colnames(summary_table))
  if (length(missing_cols) > 0) {
    stop("The following objective columns are missing from summary_table: ",
         paste(missing_cols, collapse = ", "))
  }

  if (is.null(hover_parameters)) {
    hover_parameters <- pared_result$parameter_names
  }

  missing_hover <- setdiff(hover_parameters, colnames(summary_table))
  if (length(missing_hover) > 0) {
    stop("The following hover parameter columns are missing from summary_table: ",
         paste(missing_hover, collapse = ", "))
  }

  if (is.null(xaxis_title)) xaxis_title <- x_objective
  if (is.null(yaxis_title)) yaxis_title <- y_objective

  plot_data <- data.frame(
    x = as.numeric(summary_table[[x_objective]]),
    y = as.numeric(summary_table[[y_objective]])
  )

  hover_text <- rep("", nrow(summary_table))

  for (j in seq_along(hover_parameters)) {
    cur_name <- hover_parameters[j]
    cur_val <- summary_table[[cur_name]]

    if (is.numeric(cur_val)) {
      cur_val <- round(cur_val, round_digits)
    }

    cur_piece <- paste0(cur_name, ": ", cur_val)

    if (j == 1) {
      hover_text <- cur_piece
    } else {
      hover_text <- paste0(hover_text, "<br>", cur_piece)
    }
  }

  plot_data$hover_text <- hover_text

  fig <- plotly::plot_ly(
    plot_data,
    x = ~x,
    y = ~y,
    type = "scatter",
    mode = "markers",
    name = "Pareto-optimal",
    marker = list(
      size = plot_marker_size,
      color = plot_marker_color,
      line = list(
        color = plot_marker_border_color,
        width = plot_marker_border_width
      ),
      symbol = plot_marker_symbol
    ),
    text = ~hover_text,
    hovertemplate = paste(
      paste0(xaxis_title, ": %{x}<br>"),
      paste0(yaxis_title, ": %{y}<br>"),
      "%{text}<extra></extra>"
    )
  )

  fig <- plotly::layout(
    fig,
    xaxis = list(
      title = list(text = xaxis_title, font = list(size = fontsize_x)),
      showline = show_axis_line,
      linecolor = axis_line_color,
      linewidth = axis_line_width,
      mirror = mirror_axis_line,
      showgrid = show_grid,
      gridcolor = grid_color,
      zeroline = FALSE
    ),
    yaxis = list(
      title = list(text = yaxis_title, font = list(size = fontsize_y)),
      showline = show_axis_line,
      linecolor = axis_line_color,
      linewidth = axis_line_width,
      mirror = mirror_axis_line,
      showgrid = show_grid,
      gridcolor = grid_color,
      zeroline = FALSE
    ),
    title = list(
      text = plot_title,
      font = list(size = plot_title_size),
      x = title_pos_x,
      y = title_pos_y,
      xanchor = "center",
      yanchor = "top"
    ),
    margin = list(t = 100, l = 90, r = 40, b = 90)
  )

  return(fig)
}
