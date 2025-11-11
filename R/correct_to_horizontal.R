#' Remove linear trend from Y relative to X
#'
#' This function performs a linear regression of \code{Y} on \code{X} and removes the fitted trend from \code{Y},
#' effectively flattening it to a horizontal baseline while preserving its residual variation.
#' The corrected values are returned as a new column \code{Y_corrected} in the input data frame.
#'
#' @param df A data frame containing two numeric columns: \code{X} and \code{Y}.
#'
#' @return A data frame identical to \code{df}, with an additional column \code{Y_corrected} representing
#' the detrended version of \code{Y}, centered around its original mean.
#'
#' @examples
#' set.seed(123)
#' df <- data.frame(
#'   X = 1:100,
#'   Y = 0.5 * (1:100) + rnorm(100, sd = 5)
#' )
#' df_corrected <- correct_to_horizontal(df)
#' plot(df$X, df$Y, type = "l", col = "blue")
#' lines(df$X, df_corrected$Y_corrected, col = "red")
#'
#' @export
correct_to_horizontal <- function(df) {
  set.seed(111)

  if (!all(c("X", "Y") %in% colnames(df))) {
    stop("df must include 'X' and 'Y' columns")
  }

  fit <- lm(Y ~ X, data = df)
  residuals <- resid(fit)
  y_mean <- mean(df$Y, na.rm = TRUE)
  df$Y <- residuals + y_mean

  return(df)
}
