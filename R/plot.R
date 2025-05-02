#' Plot the relationship between performance score vs. \eqn{\eta}
#'
#' This function takes a model fitted from \code{cv.coxkl} to plot the
#' relationship between performance score and corresponding \eqn{\eta}.
#'
#' @param object A model fitted from \code{cv.coxkl}.
#' @param line_size Size of the lines.
#' @param point_size Size of the points.
#' @param legend_title Title displayed above the color/linetype legend.
#' @param method_name Label used for the main method in legend, color, and linetype.
#' @param baseline_name Label used for the baseline in legend, color, and linetype.
#' @param methodline_color Color of the
#' @param baseline_color Color for the baseline point & segment
#' @param method_linetype Line type
#' @param baseline_linetype Line type
#'
#' @return
#' A ggplot object representing the relationship between performance score and eta.
#'
#'
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot geom_line geom_segment aes scale_color_manual geom_point expansion theme scale_y_continuous scale_x_continuous scale_linetype_manual element_line element_text element_blank theme_minimal
#' @method plot coxklpred
#' @export
plot.coxklpred <- function(object, line_size = 1, point_size = 3,
                    methodline_color = "#4682B4",
                    baseline_color = "#778899",
                    method_linetype  = "solid",
                    baseline_linetype= "dashed") {
  if (missing(object)) stop ("Argument 'object' is required!",call.=F)

  if (!inherits(object, "coxklpred"))
    stop("`object` must be the result of predict.coxkl().", call. = FALSE)

  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("Package \"ggplot2\" is required for this plot.", call. = FALSE)

  legend_title = ""
  method_name = "CoxKL"
  baseline_name = "Internal Only"

  col_vec <- setNames(c(methodline_color, baseline_color),
                      c(method_name, baseline_name))
  lty_vec <- setNames(c(method_linetype, baseline_linetype),
                      c(method_name, baseline_name))

  criteria <- "likelihood"
  eta_list <- object$eta_list
  performance_score <- object$likelihood

  plot_data <- data.frame(eta = eta_list, performance_score = performance_score)
  plot_data <- plot_data[order(plot_data$eta), ]

  x_value_for_dot <- plot_data[1,1]
  internal_loss1 <- plot_data[1,2]
  x_max <- max(eta_list)

  n <- nrow(plot_data)
  difference <- max(performance_score) - min(performance_score)
  loss_min <- min(performance_score)
  loss_max <- max(performance_score)

  plot_p <- ggplot() +
    geom_line(aes(x = eta, y = performance_score, color = method_name, linetype = method_name), size = line_size, data = plot_data) +
    geom_point(aes(x = x_value_for_dot, y = internal_loss1, color = baseline_name), size = point_size) +
    geom_segment(aes(x = x_value_for_dot, y = internal_loss1, xend = x_max, yend = internal_loss1 ,
                     color = baseline_name, linetype = baseline_name), size = line_size) +
    scale_y_continuous(name = criteria, limits=c(loss_min, loss_max)) +
    scale_x_continuous(name = expression(eta),  breaks = seq(min(eta_list), max(eta_list))) +
    scale_color_manual(name = legend_title, values = col_vec) +
    scale_linetype_manual(name = legend_title, values = lty_vec) +
    theme_minimal(base_size = 14) +
    theme(
      axis.line = element_line(size = 0.3),
      plot.title = element_text(size = 13),
      axis.title = element_text(size = 11),
      axis.text = element_text(size=11),
      panel.border = element_blank(),
      plot.background = element_blank(),
      panel.grid.major = element_line(color = "#d3d3d3", linewidth = 0.1),
      panel.grid.minor = element_blank()
    )

  return (plot_p)
}