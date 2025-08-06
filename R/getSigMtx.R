#' Retrieve Signature Matrix by ID
#'
#' This function extracts a subset of Signature objects from the internal Signature database
#' based on a vector of provided IDs. It validates the presence of each ID, checks for consistency
#' in detector names across selected Signatures, and returns a matrix of Signature values.
#'
#' @param ids A character vector of Signature IDs to retrieve.
#'
#' @details
#' The function reads the Signature list from the internal package path:
#' \code{system.file("sig", "Sig_list.rds", package = "USERM")}.
#' It searches for each ID in the Signature list and throws an error if any ID is not found.
#'
#' If all selected Signatures share the same set of detectors, a numeric matrix is returned
#' where rows represent detectors and columns represent Signature IDs. If detector sets differ,
#' a warning is issued and the original Signature list is returned instead.
#'
#' @return A numeric matrix of Signature values if detectors are consistent across selected Signatures.
#' Otherwise, the function returns a list of Signature objects.
#'
#' @examples
#' \dontrun{
#'   ids <- c("SCC_Bead_CD4_NFR700", "SCC_Bead_CD3_BV510",  "SCC_Bead_CD2_FITC")
#'   sig_matrix <- getSigMtx(ids)
#'   print(sig_matrix)
#' }
#'
#' @export

getSigMtx = function(ids){

  Sig_list = readRDS(system.file("sig", "Sig_list.rds", package = "USERM"))

  # find ids location
  indices = c()
  for (j in 1:length(ids)) {
    tmp_loc = NA
    for (i in 1:length(Sig_list)) {
      # print(i)
      if(Sig_list[[i]]$id == ids[j]){
        tmp_loc = i
        break
      }
    }
    if (! is.na(tmp_loc)){
      indices = c(indices,tmp_loc)
    } else {
      stop(paste0("The id: ",ids[j]," is not found in the Signature database."))
    }
  }
  Sig_list = Sig_list[indices]

  # check detectors
  first_detectors = names(Sig_list[[1]]$Signature)
  for (i in 1:length(Sig_list)) {
    tmp_detectors = names(Sig_list[[i]]$Signature)
    if(!(all((tmp_detectors %in% first_detectors)) && (length(first_detectors) == length(tmp_detectors)))){
      warning("Selected signatures have distinct detectors. Signature list is returned instead. Please check detectors and run again.")
      return(Sig_list)
    }
  }

  #prepare matrix
  mtx = matrix(nrow = length(first_detectors), ncol = length(Sig_list))
  rownames(mtx) = first_detectors
  cols = c()
  for (i in 1:ncol(mtx)) {
    mtx[(names(Sig_list[[i]]$Signature)),i] = Sig_list[[i]]$Signature
    cols = c(cols, Sig_list[[i]]$id)
  }
  colnames(mtx) = cols

  return(mtx)
}
