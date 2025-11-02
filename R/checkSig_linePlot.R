#' Plot Signature Signal Line for a Given ID
#'
#' This function generates a line plot of the signal signature for a specified ID
#' using precomputed data stored in the package. The plot shows normalized signal values
#' across detectors, with both line and point markers.
#'
#' @param id A character string specifying the ID to be visualized. Must exist in the querySig() return.
#'
#' @return A `ggplot` object showing the signal signature line plot for the given ID.
#'
#' @details The function reads the `Sig_list.rds` file from the `sig` directory of the `USERM` package.
#' If the specified ID is not found in the list, an error is raised. The plot includes a title
#' based on the ID and displays signal values across detectors.
#'
#' @examples
#' \dontrun{
#'   checkSig_linePlot("SCC_Bead_CD4_NFR700")
#' }
#'
#' @import ggplot2
#' @export


checkSig_linePlot = function(id, Sig_list = NULL){

  #check if id is in Sig_list
  if(is.null(Sig_list)){
    Sig_list = readRDS(system.file("sig", "Sig_list.rds", package = "USERM"))
  }else{
    Sig_list = Sig_list
  }

  Sig_list = readRDS(system.file("sig", "Sig_list.rds", package = "USERM"))
  if (!(id %in% names(Sig_list))) {
    stop(paste0("The input id: ", id, " is not found. Please check with querySig() first."))
  }

  Sig = Sig_list[[id]]

  df <- data.frame(Detector = factor(names(Sig$Signature),levels = names(Sig$Signature)),
                   Signal = as.numeric(Sig$Signature))
  p = ggplot(df, aes(x = Detector, y = Signal, group = 1)) +
    geom_line(color = "steelblue", size = 1) +
    geom_point(color = "darkred", size = 2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(title = Sig$id,
         x = "Detector",
         y = "Signal")

  return(p)
}
