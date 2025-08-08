#' Generate a prediction plot of unmixed spread using isobands or contour lines
#'
#' This function creates a ggplot-based visualization of predicted unmixed spread using either pseudo-color isobands or contour lines.
#' It supports multiple axis scaling methods (Linear, Log10, Arcsinh), and optionally overlays population labels based on an intensity matrix.
#'
#' @param grid A data frame containing columns \code{x}, \code{y}, and \code{z}, representing the grid coordinates and predicted values.
#' @param x_scale A character string specifying the x-axis scale. Must be one of \code{"Linear"}, \code{"Log10"}, or \code{"Arcsinh"}.
#' @param y_scale A character string specifying the y-axis scale. Must be one of \code{"Linear"}, \code{"Log10"}, or \code{"Arcsinh"}.
#' @param x_cofactor A numeric value used for Arcsinh transformation on the x-axis. Ignored for other scales.
#' @param y_cofactor A numeric value used for Arcsinh transformation on the y-axis. Ignored for other scales.
#' @param x_label A character string for the x-axis label.
#' @param y_label A character string for the y-axis label.
#' @param x_min Minimum value for the x-axis.
#' @param x_max Maximum value for the x-axis.
#' @param y_min Minimum value for the y-axis.
#' @param y_max Maximum value for the y-axis.
#' @param mode A character string specifying the plot mode. Must be either \code{"Pseudo-color"} or \code{"Contour line"}.
#' @param legned_show Logical; whether to show the legend. Default is \code{TRUE}.
#' @param label_population Optional character vector specifying population names to label on the plot.
#' @param intensity_matrix Optional matrix or data frame containing intensity values for labeling populations. Must include rows named by \code{x_label} and \code{y_label}.
#'
#' @return A \code{ggplot} object representing the prediction plot.
#'
#' @examples
#' \dontrun{
#' grid <- expand.grid(x = seq(0, 1000, length.out = 50),
#'                     y = seq(0, 1000, length.out = 50))
#' grid$z <- with(grid, exp(-((x - 500)^2 + (y - 500)^2) / 1e5))
#' PredIsobandPlot(grid, x_scale = "Linear", y_scale = "Linear",
#'                 x_cofactor = 5, y_cofactor = 5,
#'                 x_label = "X", y_label = "Y",
#'                 x_min = 0, x_max = 1000,
#'                 y_min = 0, y_max = 1000,
#'                 mode = "Pseudo-color")
#' }
#' @export
#' @importFrom ggplot2 ggplot aes geom_contour geom_contour_filled geom_text
#' @importFrom ggplot2 scale_color_viridis_c scale_x_continuous scale_y_continuous
#' @importFrom ggplot2 labs theme_light theme xlim ylim after_stat
#' @importFrom scales trans_new



PredIsobandPlot = function(grid,x_scale,y_scale,x_cofactor,y_cofactor,
                           x_label,y_label,x_min,x_max,y_min,y_max,
                           mode,legned_show = TRUE,label_population = NA,intensity_matrix = NA){
  asinh_trans_cofactor <- function(cofactor = 5) {
    scales::trans_new(
      name = paste0("asinh_", cofactor),
      transform = function(x) asinh(x / cofactor),
      inverse = function(x) cofactor * sinh(x),
      domain = NULL
    )
  }
  signed_log10_trans <- function() {
    scales::trans_new(
      name = "signed_log10",
      transform = function(x) {
        sign(x) * log10(abs(x) + 1)  # +1 to avoid log(0)
      },
      inverse = function(x) {
        sign(x) * (10^abs(x) - 1)
      },
      domain = c(-Inf, Inf)
    )
  }

  p <- ggplot(grid, aes(x = x, y = y, z = z))

  if(mode == "Pseudo-color"){
    p <- p + geom_contour_filled()
  }else if(mode == "Contour line"){
    p <- p + geom_contour(aes(color = after_stat(level)))+
      scale_color_viridis_c(option = "C")
  }

  p <- p + labs(title = "Prediction of Unmixed Spread",
         x = x_label,
         y = y_label,
         fill = "Density")
  if (x_scale == "Linear"){
    p <- p + xlim(x_min, x_max)
  }
  if (y_scale == "Linear"){
    p <- p + ylim(y_min, y_max)
  }


  # scale
  if (x_scale == "Log10") {
    p <- p + scale_x_continuous(trans = signed_log10_trans(),
                                limits = c(x_min, x_max))
  } else if (x_scale == "Arcsinh") {
    p <- p + scale_x_continuous(trans = asinh_trans_cofactor(x_cofactor),
                                limits = c(x_min, x_max))
  }

  if (y_scale == "Log10") {
    p <- p + scale_y_continuous(trans = signed_log10_trans(),
                                limits = c(y_min, y_max))
  } else if (y_scale == "Arcsinh") {
    p <- p + scale_y_continuous(trans = asinh_trans_cofactor(y_cofactor),
                                limits = c(y_min, y_max))
  }


  # add population label
  if (!is.na(label_population[1])) {
    label_df <- intensity_matrix[c(x_label,y_label),label_population]
    rownames(label_df) = c("x","y")
    label_df = as.data.frame(t(label_df))
    label_df$population = rownames(label_df)
    label_df$z = 0
    p <- p + geom_text(data = label_df, aes(x = x, y = y, label = population), vjust = -0.5)
  }


  p <- p + theme_light()

  if(!legned_show){
    p <- p + theme(legend.position = "none")
  }
  return(p)

}

