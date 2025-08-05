#' Interpolate a color based on a numeric value
#'
#' This function interpolates a color between three specified colors (minimum, midpoint, and maximum)
#' based on the position of a numeric value within a defined range.
#'
#' @param val A numeric value to be mapped to a color.
#' @param min The minimum value of the range.
#' @param mid The midpoint value of the range.
#' @param max The maximum value of the range.
#' @param colormin The color (in HEX format) corresponding to the minimum value.
#' @param colormid The color (in HEX format) corresponding to the midpoint value.
#' @param colormax The color (in HEX format) corresponding to the maximum value.
#'
#' @return A HEX color string representing the interpolated color.
#' @examples
#' getColor(0.75, 0, 0.5, 1, "#0000FF", "#00FF00", "#FF0000")
#' @export

getColor <- function(val, min, mid, max, colormin, colormid, colormax) {
  # Clamp val within [min, max]
  val <- max(min, min(max, val))

  # Convert HEX to RGB vector
  hex_to_rgb <- function(hex) {
    rgb <- col2rgb(hex) / 255
    return(as.numeric(rgb))
  }

  rgb_min <- hex_to_rgb(colormin)
  rgb_mid <- hex_to_rgb(colormid)
  rgb_max <- hex_to_rgb(colormax)

  # Interpolate
  if (val <= mid) {
    ratio <- (val - min) / (mid - min)
    rgb <- rgb_min + ratio * (rgb_mid - rgb_min)
  } else {
    ratio <- (val - mid) / (max - mid)
    rgb <- rgb_mid + ratio * (rgb_max - rgb_mid)
  }

  # Convert RGB to HEX
  rgb_to_hex <- function(rgb) {
    rgb(rgb[1], rgb[2], rgb[3])
  }

  return(rgb_to_hex(rgb))
}
