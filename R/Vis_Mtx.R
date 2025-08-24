#' Visualize a Matrix as a Heatmap with Custom Colors and Labels
#'
#' This function creates a heatmap visualization of a numeric matrix using the
#' `ComplexHeatmap` package. It allows customization of the color scale and
#' overlays numeric values on each cell.
#'
#' @param mat A numeric matrix to be visualized. Must not be `NULL`.
#' @param min Minimum value for the color scale. Default is -1.
#' @param mid Midpoint value for the color scale. Default is 0.
#' @param max Maximum value for the color scale. Default is 1.
#' @param mincolor Color representing the minimum value. Default is `"blue"`.
#' @param midcolor Color representing the midpoint value. Default is `"white"`.
#' @param maxcolor Color representing the maximum value. Default is `"red"`.
#'
#' @return A `Heatmap` object representing the input matrix with customized colors
#' and cell labels.
#'
#' @details The function uses `colorRamp2` to define a color gradient and
#' `ComplexHeatmap::Heatmap` to generate the heatmap. Each cell displays its
#' numeric value formatted to one decimal place. Row and column clustering are disabled.
#'
#' @examples
#' \dontrun{
#' mat <- matrix(runif(25, -1, 1), nrow = 5)
#' Vis_Mtx(mat)
#' }
#'
#' @import ComplexHeatmap
#' @import circlize
#' @import grid
#' @export


Vis_Mtx <- function(mat = NULL, min = -1, mid = 0, max = 1, mincolor = "blue", midcolor = "white", maxcolor = "red") {
  if (is.null(mat) || !is.matrix(mat)) {
    stop("mat is missing or not a matrix.")
  }

  col_fun = colorRamp2(c(min, mid, max), c(mincolor, midcolor, maxcolor))

  p = Heatmap(mat,col = col_fun, cluster_rows = FALSE,cluster_columns = FALSE,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 6))
              },column_title = paste0("slopMtx of ", Res$id),name = "beta")

  return(p)
}
