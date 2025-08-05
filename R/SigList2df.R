#' Convert a list of signature objects into a data frame
#'
#' This function takes a list of signature objects (each containing a named numeric vector under \code{Signature})
#' and combines them into a single data frame. Each column in the resulting data frame represents one signature,
#' and the column names are assigned based on a specified identifier field.
#'
#' @param SigList A list of signature objects, each containing an \code{id}, \code{PrimaryName}, \code{SecondaryName}, and a named numeric vector \code{Signature}.
#' @param SigName A character string specifying which identifier to use as column names in the output data frame.
#'        Must be one of \code{"id"}, \code{"PrimaryName"}, or \code{"SecondaryName"}. Default is \code{"id"}. "id" is recommended for consistence in the following use.
#'
#' @return A data frame where each column corresponds to a signature, and each row corresponds to a detector channel.
#'
#' @examples
#' \dontrun{
#' sig1 <- list(id = "FITC_CD3", PrimaryName = "FITC", SecondaryName = "CD3",
#'              Signature = c(D1 = 0.2, D2 = 0.8))
#' sig2 <- list(id = "PE_CD4", PrimaryName = "PE", SecondaryName = "CD4",
#'              Signature = c(D1 = 0.5, D2 = 0.5))
#' sig_df <- SigList2df(list(sig1, sig2), SigName = "id")
#' }
#' @export

SigList2df = function(SigList,SigName = "id"){
  for (i in 1:length(SigList)) {
    if(SigName == "id"){
      sigid_tmp = SigList[[i]]$id
    }else if(SigName == "PrimaryName"){
      sigid_tmp = SigList[[i]]$PrimaryName
    }else if(SigName == "SecondaryName"){
      sigid_tmp = SigList[[i]]$SecondaryName
    }else{
      stop(paste0("SigName needs to be 'id', 'PrimaryName', or 'SecondaryName'"))
    }

    if(i == 1){
      cols = names(SigList[[1]]$Signature)
      SigDf = as.data.frame(SigList[[i]]$Signature)
      colnames(SigDf)[i] = sigid_tmp
    }else{
      #check detector names
      tmp_cols = names(SigList[[i]]$Signature)
      unmatched <- all(cols == tmp_cols)
      if (!unmatched) {
        stop(paste("Error: The following Signatures do not have the same detectors:",
                   paste(c(names(SigList)[1],names(SigList)[i]), collapse = ", ")))
      }
      #merge sig
      SigDf = cbind(SigDf,as.data.frame(SigList[[i]]$Signature))
      colnames(SigDf)[i] = sigid_tmp
    }
  }
  return(SigDf)
}
