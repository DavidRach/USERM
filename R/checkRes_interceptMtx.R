
#' Visualize the Intercept Matrix from a Result Object
#'
#' This function creates a heatmap visualization of the `interceptMtx` matrix contained within a result object `Res`.
#' It uses a customizable color gradient and overlays numeric values on each cell.
#'
#' @param Res A Res object.
#' @param min Numeric value representing the minimum of the color scale. Default is -10.
#' @param mid Numeric value representing the midpoint of the color scale. Default is 0.
#' @param max Numeric value representing the maximum of the color scale. Default is 10.
#' @param mincolor Color corresponding to the minimum value. Default is "blue".
#' @param midcolor Color corresponding to the midpoint value. Default is "white".
#' @param maxcolor Color corresponding to the maximum value. Default is "red".
#'
#' @return A `Heatmap` object from the ComplexHeatmap package, visualizing the `interceptMtx` with annotated values.
#'
#' @examples
#' \dontrun{
#'   checkRes_interceptMtx(result)
#' }
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#' @export

checkRes_interceptMtx <- function(Res, min = -10, mid = 0, max = 10, mincolor = "blue", midcolor = "white", maxcolor = "red") {
  if (is.null(Res$interceptMtx) || !is.matrix(Res$interceptMtx)) {
    stop("Res$interceptMtx is missing or not a matrix.")
  }
  mat = Res$interceptMtx
  col_fun = colorRamp2(c(min, mid, max), c(mincolor, midcolor, maxcolor))

  p = Heatmap(mat,col = col_fun, cluster_rows = FALSE,cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
              },column_title = paste0("interceptMtx of ", Res$id),name = "intercept")

  return(p)
}
