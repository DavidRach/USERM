
#' Plot Covariance Scatter Between Two Detectors with Regression Line
#'
#' This function generates a scatter plot of covariance values between two specified detectors
#' across bin midpoints, using data from a result object `Res`. A regression line is added
#' using the corresponding slope and intercept values, which are also displayed in the plot title.
#'
#' @param Res A list containing the following components:
#' \itemize{
#'   \item \code{bin_mids}: A numeric vector of bin midpoints.
#'   \item \code{detectors}: A character vector of detector names.
#'   \item \code{cov_matrices}: A 3D array of covariance values indexed by detector pairs and bins.
#'   \item \code{slopMtx}: A matrix of slope values for each detector pair.
#'   \item \code{interceptMtx}: A matrix of intercept values for each detector pair.
#' }
#' @param detector1 Name of the one detector (character string).
#' @param detector2 Name of the another detector (character string).
#'
#' @return A `ggplot` object showing the scatter plot of covariance values and a regression line.
#'
#' @examples
#' \dontrun{
#' checkRes_covScatter(Res = ResObj,
#' detector1 = ResObj$detectors[1],
#' detector2 = ResObj$detectors[2])
#' }

#' @import ggplot2
#' @export
checkRes_covScatter <- function(Res, detector1, detector2) {

  if (is.null(Res$cov_matrices) || !is.matrix(Res$cov_matrices)) {
    stop("Res$cov_matrices is missing or not a matrix.")
  }

  x <- Res$bin_mids
  detectors <- Res$detectors

  # Find the index of the specified detectors
  idx1 <- which(detectors == detector1)
  idx2 <- which(detectors == detector2)

  if (length(idx1) == 0 || length(idx2) == 0) {
    stop("Specified detectors not found in Res$detectors.")
  }

  y <- Res$cov_matrices[idx1, idx2, ]
  slope <- Res$slopMtx[idx1, idx2]
  intercept <- Res$interceptMtx[idx1, idx2]

  df <- data.frame(x = x, y = y)

  title_text <- sprintf(
    "Covariance between %s and %s (slope = %.3f, intercept = %.3f)",
    detector1, detector2, slope, intercept
  )

  p <- ggplot(df, aes(x = x, y = y)) +
    geom_point(color = "black") +
    geom_abline(slope = slope, intercept = intercept, color = "red", linetype = "dashed") +
    labs(
      title = title_text,
      x = "Bin Midpoints",
      y = "Covariance"
    ) +
    theme_light()

  return(p)
}
