#' Add a Res object to a Userm object
#'
#' This function validates and adds a Res object (`Res`) to a Userm object (`Userm`). It checks whether the `Res` object has a valid ID, matching detectors, and properly structured matrices (`interceptMtx` and `slopMtx`). If all checks pass, the Res object is added to the `Userm$Res` list.
#'
#' @param Res A Res object.
#' @param Userm A Userm object, generated with CreateUserm function.
#'
#' @return The updated `Userm` object with the new `Res` object added.
#' @examples
#' \dontrun{
#' Res <- CreateRes(id, R, A)
#' Userm <- CreateUserm(A)
#' Userm <- AddRes2Userm(Res, Userm)
#' }
#' @export

AddRes2Userm = function(Res, Userm){

  fluors = Userm$fluors
  detectors = Userm$detectors

  #check if id of Res is in fluors of Userm
  Res_id = Res$id
  if(!(Res_id %in% fluors)){
    stop(paste0("Error: the input Res (id: ",Res_id,") is not in the Userm."))
  }

  #check Res$detectors
  if(!all(detectors == Res$detectors)){
    stop(paste0("Error: the input Res has inconsistent detectors with that in the Userm."))
  }

  #check Res$interceptMtx

  if(is.null(Res$interceptMtx) || is.na(Res$interceptMtx)[1]){
    stop(paste0("Error: the input Res has no interceptMtx. Please do SlopEstimation first."))
  }
  if(!all(dim(Res$interceptMtx) == c(length(detectors), length(detectors)))){
    stop(paste0("Error: the interceptMtx of input Res is not a (",length(detectors),", ",length(detectors),") matrix.  Please check it first."))
  }

  #check Res$slopMtx
  if(is.null(Res$slopMtx) || is.na(Res$slopMtx)[1]){
    stop(paste0("Error: the input Res has no slopMtx Please do SlopEstimation first."))
  }
  if(!all(dim(Res$slopMtx) == c(length(detectors), length(detectors)))){
    stop(paste0("Error: the slopMtx of input Res is not a (",length(detectors),", ",length(detectors),") matrix.  Please check it first."))
  }

  #if all good
  Userm$Res[[Res_id]] = Res

  return(Userm)
}
