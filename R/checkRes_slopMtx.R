#' Visualize the Slope Matrix from a Result Object
#'
#' This function generates a heatmap visualization of the `slopMtx` matrix contained within a result object `Res`.
#' It uses a customizable color gradient and overlays numeric values on each cell.
#'
#' @param Res A Res object.
#' @param min Numeric value representing the minimum of the color scale. Default is -1.
#' @param mid Numeric value representing the midpoint of the color scale. Default is 0.
#' @param max Numeric value representing the maximum of the color scale. Default is 1.
#' @param mincolor Color corresponding to the minimum value. Default is "blue".
#' @param midcolor Color corresponding to the midpoint value. Default is "white".
#' @param maxcolor Color corresponding to the maximum value. Default is "red".
#'
#' @return A `Heatmap` object from the ComplexHeatmap package, visualizing the `slopMtx` with annotated values.
#'
#' @examples
#' \dontrun{
#'   checkRes_slopMtx(Res)
#' }
#'
#' @import ComplexHeatmap
#' @import circlize
#' @export

checkRes_slopMtx <- function(Res, min = -1, mid = 0, max = 1, mincolor = "blue", midcolor = "white", maxcolor = "red") {
  if (is.null(Res$slopMtx) || !is.matrix(Res$slopMtx)) {
    stop("Res$slopMtx is missing or not a matrix.")
  }
  mat = Res$slopMtx
  col_fun = colorRamp2(c(min, mid, max), c(mincolor, midcolor, maxcolor))

  p = Heatmap(mat,col = col_fun, cluster_rows = FALSE,cluster_columns = FALSE,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
          },column_title = paste0("slopMtx of ", Res$id),name = "beta")

  return(p)
}
