#' Overlay Prediction Grids Using Isobands or Contour Lines
#'
#' This function creates a ggplot-based visualization that overlays multiple prediction grids
#' on a single plot. Each grid represents a model's prediction and is distinguished by color.
#' The function supports various axis scaling methods (Linear, Log10, Arcsinh), and optionally
#' overlays population labels based on an intensity matrix.
#'
#' @param grid_list A named list of data frames, each containing columns \code{x}, \code{y}, and \code{z}.
#'                  Each data frame represents a prediction grid for a model.
#' @param x_label A character string for the x-axis label.
#' @param y_label A character string for the y-axis label.
#' @param x_scale A character string specifying the x-axis scale. Must be one of \code{"Linear"}, \code{"Log10"}, or \code{"Arcsinh"}.
#' @param y_scale A character string specifying the y-axis scale. Must be one of \code{"Linear"}, \code{"Log10"}, or \code{"Arcsinh"}.
#' @param x_cofactor A numeric value used for Arcsinh transformation on the x-axis. Ignored for other scales.
#' @param y_cofactor A numeric value used for Arcsinh transformation on the y-axis. Ignored for other scales.
#' @param x_min Minimum value for the x-axis.
#' @param x_max Maximum value for the x-axis.
#' @param y_min Minimum value for the y-axis.
#' @param y_max Maximum value for the y-axis.
#' @param mode A character string specifying the plot mode. Must be either \code{"Contour line"} or \code{"Pseudo-color"}.
#' @param legend_show Logical; whether to show the legend. Default is \code{TRUE}.
#' @param label_population Optional character vector specifying population names to label on the plot.
#' @param intensity_matrix Optional matrix or data frame containing intensity values for labeling populations.
#'                         Must include rows named by \code{x_label} and \code{y_label}.
#'
#' @return A \code{ggplot} object representing the combined prediction plot.
#'
#' @examples
#' \dontrun{
#' PredIsobandPlotMulti(grid_list,
#'                      x_label = "X", y_label = "Y",
#'                      x_scale = "Linear", y_scale = "Linear",
#'                      x_min = 0, x_max = 1000,
#'                      y_min = 0, y_max = 1000,
#'                      mode = "Pseudo-color")
#' }
#' @export
#' @import ggplot2
#' @import scales
#' @import dplyr


PredIsobandPlotMulti <- function(grid_list,
                                 x_label = "X", y_label = "Y",
                                 x_scale = "Linear", y_scale = "Linear",
                                 x_cofactor = 5, y_cofactor = 5,
                                 x_min = NULL, x_max = NULL,
                                 y_min = NULL, y_max = NULL,
                                 mode = "Pseudo-color",
                                 legend_show = TRUE,
                                 label_population = NA,
                                 intensity_matrix = NA) {

  library(ggplot2)
  library(scales)
  library(dplyr)

  # Custom transformation functions
  asinh_trans_cofactor <- function(cofactor = 5) {
    trans_new(
      name = paste0("asinh_", cofactor),
      transform = function(x) asinh(x / cofactor),
      inverse = function(x) cofactor * sinh(x)
    )
  }

  signed_log10_trans <- function() {
    trans_new(
      name = "signed_log10",
      transform = function(x) sign(x) * log10(abs(x) + 1),
      inverse = function(x) sign(x) * (10^abs(x) - 1)
    )
  }

  # Combine all grids and add model name
  combined_grid <- bind_rows(
    lapply(names(grid_list), function(name) {
      df <- grid_list[[name]]
      df$model <- name
      return(df)
    })
  )

  # Initialize ggplot
  p <- ggplot(combined_grid, aes(x = x, y = y, z = z))

  # Add plot layers
  if (mode == "Pseudo-color") {
    p <- p + geom_contour_filled(aes(fill = model), alpha = 0.1)
  } else if (mode == "Contour line") {
    p <- p + geom_contour(aes(color = model))
  }

  # Axis labels
  p <- p + labs(title = "Overlay of Multiple Prediction Grids",
                x = x_label, y = y_label)

  # Axis scaling
  if (!is.null(x_min) && !is.null(x_max)) {
    if (x_scale == "Linear") {
      p <- p + xlim(x_min, x_max)
    } else if (x_scale == "Log10") {
      p <- p + scale_x_continuous(trans = signed_log10_trans(), limits = c(x_min, x_max))
    } else if (x_scale == "Arcsinh") {
      p <- p + scale_x_continuous(trans = asinh_trans_cofactor(x_cofactor), limits = c(x_min, x_max))
    }
  }

  if (!is.null(y_min) && !is.null(y_max)) {
    if (y_scale == "Linear") {
      p <- p + ylim(y_min, y_max)
    } else if (y_scale == "Log10") {
      p <- p + scale_y_continuous(trans = signed_log10_trans(), limits = c(y_min, y_max))
    } else if (y_scale == "Arcsinh") {
      p <- p + scale_y_continuous(trans = asinh_trans_cofactor(y_cofactor), limits = c(y_min, y_max))
    }
  }

  # Add population label
  if (!is.na(label_population[1])) {
    label_df <- intensity_matrix[c(x_label, y_label), label_population]
    rownames(label_df) <- c("x", "y")
    label_df <- as.data.frame(t(label_df))
    label_df$population <- rownames(label_df)
    label_df$z <- 0
    p <- p + geom_text(data = label_df, aes(x = x, y = y, label = population), vjust = -0.5)
  }

  # Theme and legend
  p <- p + theme_light()
  if (!legend_show) {
    p <- p + theme(legend.position = "none")
  }

  return(p)
}
