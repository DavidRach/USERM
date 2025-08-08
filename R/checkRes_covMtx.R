#' Visualize Covariance Matrix at a Specific Bin
#'
#' This function generates a heatmap of the covariance matrix at a specified bin index
#' from a result object `Res`. The heatmap uses a customizable color gradient and overlays
#' numeric values on each cell.
#'
#' @param Res A Res object containing the following components:
#' \itemize{
#'   \item \code{cov_matrices}: A 3D array of covariance matrices (detector x detector x bin).
#'   \item \code{detectors}: A character vector of detector names.
#'   \item \code{bin_mids}: A numeric vector of bin midpoints.
#'   \item \code{id}: An identifier for the result object, used in the plot title.
#' }
#' @param bin Integer index specifying which bin's covariance matrix to visualize. Default is 1 (the first bin).
#' @param min Minimum value for the color scale. Default is -1.
#' @param mid Midpoint value for the color scale. Default is 0.
#' @param max Maximum value for the color scale. Default is 1.
#' @param mincolor Color corresponding to the minimum value. Default is "blue".
#' @param midcolor Color corresponding to the midpoint value. Default is "white".
#' @param maxcolor Color corresponding to the maximum value. Default is "red".
#'
#' @return A `Heatmap` object from the ComplexHeatmap package, visualizing the covariance matrix.
#'
#' @examples
#' \dontrun{
#'   checkRes_covMtx(Res, bin = 2)
#' }
#'
#' @import ComplexHeatmap
#' @import circlize
#' @export


checkRes_covMtx <- function(Res, bin = 1, min = -1, mid = 0, max = 1, mincolor = "blue", midcolor = "white", maxcolor = "red") {
  if (is.null(Res$cov_matrices[,,bin]) || !is.matrix(Res$cov_matrices[,,bin])) {
    stop("Res$cov_matrices[,,bin] is missing or not a matrix.")
  }
  mat = Res$cov_matrices[,,bin]
  rownames(mat) = Res$detectors
  colnames(mat) = Res$detectors

  col_fun = colorRamp2(c(min, mid, max), c(mincolor, midcolor, maxcolor))

  p = Heatmap(mat,col = col_fun, cluster_rows = FALSE,cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
              },column_title = paste0("cov_matrix of ", Res$id, " at bin_mid = ",round(Res$bin_mids[bin],digits = 2)),name = "cov")

  return(p)
}
