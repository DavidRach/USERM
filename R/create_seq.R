#' Generate a numeric sequence based on a scaling method
#'
#' This function generates a numeric sequence between a minimum and maximum value using one of three scaling methods:
#' \code{"Linear"}, \code{"Log10"}, or \code{"Arcsinh"}. It supports signed transformations for logarithmic and arcsinh scaling.
#'
#' @param min A numeric value specifying the start of the sequence.
#' @param max A numeric value specifying the end of the sequence.
#' @param len An integer specifying the number of points in the sequence.
#' @param scale A character string specifying the scaling method. Must be one of \code{"Linear"}, \code{"Log10"}, or \code{"Arcsinh"}.
#' @param cofactor A numeric value used in the \code{"Arcsinh"} transformation. Ignored for other scaling methods.
#'
#' @return A numeric vector of length \code{len} representing the scaled sequence.
#'
#' @examples
#' create_seq(min = 0, max = 1000, len = 50, scale = "Linear", cofactor = 5)
#' create_seq(min = -100, max = 1000, len = 50, scale = "Log10", cofactor = 5)
#' create_seq(min = -100, max = 1000, len = 50, scale = "Arcsinh", cofactor = 5)
#'
#' @export


create_seq <- function(min, max, len, scale, cofactor) {
  if (scale == "Linear") {
    output <- seq(min, max, length.out = len)

  } else if (scale == "Log10") {
    # use signed log10 transform
    min_trans <- sign(min) * log10(abs(min) + 1)
    max_trans <- sign(max) * log10(abs(max) + 1)
    seq_trans <- seq(min_trans, max_trans, length.out = len)
    output <- sign(seq_trans) * (10^abs(seq_trans) - 1)

  } else if (scale == "Arcsinh") {
    min_trans <- asinh(min / cofactor)
    max_trans <- asinh(max / cofactor)
    seq_trans <- seq(min_trans, max_trans, length.out = len)
    output <- sinh(seq_trans) * cofactor

  } else {
    stop("Unsupported scale type") #will not be used
  }

  return(output)
}
