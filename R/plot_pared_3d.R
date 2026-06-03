#' Three-Dimensional Pareto Front Plot
#'
#' Creates a three-dimensional Plotly scatter plot from a \code{pared_result}
#' object returned by \code{\link{pared_optimize}}. The user selects three
#' objective columns from the Pareto summary table to display on the x-, y-, and
#' z-axes. Optional projection lines can be drawn from each Pareto point to the
#' xy-plane.
#'
#' @param pared_result An object of class \code{"pared_result"}, typically
#'   returned by \code{\link{pared_optimize}}.
#' @param x_objective Character string giving the name of the objective column
#'   to display on the x-axis.
#' @param y_objective Character string giving the name of the objective column
#'   to display on the y-axis.
#' @param z_objective Character string giving the name of the objective column
#'   to display on the z-axis.
#' @param hover_parameters Character vector giving the parameter columns to
#'   display in the hover text. Default is \code{NULL}, in which case
#'   \code{pared_result$parameter_names} is used.
#' @param plot_title Character string giving the plot title. Default is
#'   \code{"Pareto front"}.
#' @param xaxis_title Character string giving the x-axis title. Default is
#'   \code{NULL}, in which case \code{x_objective} is used.
#' @param yaxis_title Character string giving the y-axis title. Default is
#'   \code{NULL}, in which case \code{y_objective} is used.
#' @param zaxis_title Character string giving the z-axis title. Default is
#'   \code{NULL}, in which case \code{z_objective} is used.
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
#' @param plot_title_size Font size for the plot title. Default is \code{30}.
#' @param title_pos_x Horizontal plot-title position on the Plotly 0--1 scale.
#'   Default is \code{0.5}.
#' @param title_pos_y Vertical plot-title position on the Plotly 0--1 scale.
#'   Default is \code{0.90}.
#' @param draw_projection Numeric flag, either \code{0} or \code{1}, indicating
#'   whether projection lines from each point to the xy-plane are drawn.
#'   Default is \code{1}.
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
#' @param round_digits Number of digits used when rounding hover-label values.
#'   Default is \code{3}.
#'
#' @return A \code{plotly} object containing a three-dimensional Pareto-front
#'   scatter plot.
#'
#' @details
#' The function extracts the selected objective columns from
#' \code{pared_result$summary_table}. Hover text is constructed from the selected
#' parameter columns. When \code{draw_projection = 1}, dashed projection lines
#' are drawn from each Pareto point to the minimum observed z-plane. This function
#' is useful for visualizing three-objective Pareto fronts or three-dimensional
#' projections of higher-dimensional Pareto results.
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
#' toy_res$summary_table
#'
#' plot_pared_3d(toy_res,
#'               x_objective = "Mean",
#'               y_objective = "Harmonic mean",
#'               z_objective = "MAD around 1.5",
#'               plot_title = "Toy Pareto front",
#'               plot_marker_color = "blue"
#'               )
#'               
#'}
#'
#' @seealso \code{\link{pared_optimize}}, \code{\link{plot_pared_2d}}
#'
#' @export


plot_pared_3d <- function(pared_result,
                          x_objective,
                          y_objective,
                          z_objective,
                          hover_parameters = NULL,
                          plot_title = "Pareto front",
                          xaxis_title = NULL,
                          yaxis_title = NULL,
                          zaxis_title = NULL,
                          plot_marker_color = 'rgba(255, 0, 0, 0.7)',
                          plot_marker_border_color = 'black',
                          plot_marker_size = 10,
                          plot_marker_border_width = 2,
                          plot_marker_symbol = 'circle',
                          fontsize_x = 20,
                          fontsize_y = 20,
                          fontsize_z = 20,
                          plot_title_size = 30,
                          title_pos_x = 0.5,
                          title_pos_y = 0.90,
                          draw_projection = 1,
                          proj_line_color = 'black',
                          proj_line_width = 2,
                          proj_marker_size = 2,
                          proj_marker_color = 'rgba(255, 0, 0, 0.7)',
                          proj_marker_border_width = 2,
                          proj_marker_border_color = 'black',
                          round_digits = 3) {
  
  if (!inherits(pared_result, "pared_result")) {
    stop("pared_result must be an object returned by pared_optimize().")
  }
  
  summary_table <- pared_result$summary_table
  
  required_cols <- c(x_objective, y_objective, z_objective)
  
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
  if (is.null(zaxis_title)) zaxis_title <- z_objective
  
  plot_data <- data.frame(
    x = as.numeric(summary_table[[x_objective]]),
    y = as.numeric(summary_table[[y_objective]]),
    z = as.numeric(summary_table[[z_objective]])
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
    z = ~z,
    type = 'scatter3d',
    mode = 'markers',
    name = 'Pareto-optimal',
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
      paste0(zaxis_title, ": %{z}<br>"),
      "%{text}<extra></extra>"
    )
  )
  
  fig <- plotly::layout(
    fig,
    scene = list(
      xaxis = list(title = list(text = xaxis_title, font = list(size = fontsize_x))),
      yaxis = list(title = list(text = yaxis_title, font = list(size = fontsize_y))),
      zaxis = list(title = list(text = zaxis_title, font = list(size = fontsize_z)))
    ),
    title = list(
      text = plot_title,
      font = list(size = plot_title_size),
      x = title_pos_x,
      y = title_pos_y,
      xanchor = "center",
      yanchor = "top"
    ),
    margin = list(t = 50, l = 0, r = 0, b = 0)
  )
  
  if (draw_projection == 1) {
    z_min <- min(plot_data$z, na.rm = TRUE)
    
    for (i in seq_len(nrow(plot_data))) {
      fig <- plotly::add_trace(
        fig,
        x = c(plot_data$x[i], plot_data$x[i]),
        y = c(plot_data$y[i], plot_data$y[i]),
        z = c(z_min, plot_data$z[i]),
        text = plot_data$hover_text[i],
        type = 'scatter3d',
        mode = 'lines + markers',
        showlegend = FALSE,
        line = list(
          color = proj_line_color,
          width = proj_line_width,
          dash = 'dash'
        ),
        marker = list(
          size = proj_marker_size,
          color = proj_marker_color,
          line = list(
            color = proj_marker_border_color,
            width = proj_marker_border_width
          )
        )
      )
    }
  }
  
  return(fig)
}