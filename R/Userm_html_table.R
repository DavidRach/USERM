#' Generate an HTML table with color-coded cells and a legend
#'
#' This function creates an HTML representation of a numeric matrix with each cell
#' color-coded based on its value. A color legend is also appended to help interpret
#' the color scale.
#'
#' @param mat A numeric matrix to be displayed as an HTML table.
#' @param val_min Minimum value for color scaling.
#' @param val_mid Midpoint value for color scaling.
#' @param val_max Maximum value for color scaling.
#' @param colormin Color corresponding to `val_min`.
#' @param colormid Color corresponding to `val_mid`.
#' @param colormax Color corresponding to `val_max`.
#' @param Legend A character string for the legend title. Default is `"Color Legend"`.
#'
#' @return A character string containing HTML code for the table and legend.
#'
#' @details This function is useful for visualizing matrices with a heatmap-like
#' color gradient directly in Shiny apps or HTML reports. It relies on a helper
#' function `getColor()` to compute the appropriate background color for each cell.
#'
#' @examples
#' mat <- matrix(runif(9, 0, 100), nrow = 3)
#' rownames(mat) <- paste0("Row", 1:3)
#' colnames(mat) <- paste0("Col", 1:3)
#' Userm_html_table(mat, 0, 50, 100, "#0000FF", "#FFFFFF", "#FF0000")
#'
#' @export


Userm_html_table = function(mat,val_min, val_mid, val_max, colormin, colormid, colormax, Legend = "Color Legend"){
  # construct HTML
  html <- "<table style='border-collapse:collapse;'>"
  html <- paste0(html, "<tr><th></th>",
                 paste0("<th style='padding: 6px 12px;'>", colnames(mat), "</th>", collapse = ""),
                 "</tr>")

  for (i in seq_along(rownames(mat))) {
    html <- paste0(html, "<tr><th>", rownames(mat)[i], "</th>")
    for (j in seq_along(colnames(mat))) {
      val <- mat[i, j]
      color <-getColor(val = val, min = val_min, mid = val_mid, max = val_max,
                       colormin = colormin, colormid = colormid, colormax = colormax)

      html <- paste0(html, sprintf(
        "<td style='background-color:%s; text-align:center; padding:4px;'>%.2f</td>",
        color, val
      ))
    }
    html <- paste0(html, "</tr>")
  }
  html <- paste0(html, "</table>")

  legend_html <- sprintf("
                          <div style='margin-top:15px; font-family:Arial; font-size:12px;'>
                            <div style='margin-bottom:4px;'>%s</div>
                            <div style='width:50%%; display:flex; justify-content:space-between; margin-top:4px;'>
                              <span >%.2f</span>
                              <span >%.2f</span>
                              <span >%.2f</span>
                            </div>
                            <div style='width:50%%; height:20px; background:linear-gradient(to right, %s, %s, %s); border:1px solid #ccc;'></div>

                          </div>
                        ",Legend,
                         val_min, val_mid, val_max,
                         colormin, colormid, colormax)
  final_html = paste0(html, legend_html)
  return(final_html)
}
