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
