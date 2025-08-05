#' Add scaling parameters for a specific fluorophore
#'
#' This function adds or updates scaling parameters for a given fluorophore in a userm object.
#' It validates the input fluor name, scaling method, and numeric parameters, then stores them in the
#' `Scale_df` field of the `Userm` object.
#'
#' @param Userm A userm object created by \code{\link{CreateUserm}}, containing fluor and detector information.
#' @param fluor A character string specifying the fluorophore to update.
#' @param min A numeric value specifying the minimum intensity threshold. Default is -100.
#' @param max A numeric value specifying the maximum intensity threshold. Default is 1000.
#' @param scale A character string specifying the scaling method. Must be one of \code{"Linear"}, \code{"Log10"}, or \code{"Arcsinh"}. Default is \code{"Linear"}.
#' @param cofactor A positive numeric value used for scaling (only used for \code{"Arcsinh"}). Default is 10.
#'
#' @return The updated \code{Userm} object with new scaling parameters for the specified fluorophore.
#'
#' @examples
#' \dontrun{
#' A <- matrix(runif(20), nrow = 4, ncol = 5)
#' rownames(A) <- paste0("D", 1:4)
#' colnames(A) <- paste0("F", 1:5)
#' Userm <- CreateUserm(A)
#' Userm <- AddScale(Userm, fluor = "F1", min = 0, max = 5000, scale = "Arcsinh", cofactor = 5)
#' }
#' @export


AddScale = function(Userm,fluor,min = -100, max = 1000, scale = "Linear", cofactor = 10){
  # scale = c("Linear","Log10","Arcsinh")
  fluors = Userm$fluors

  #check if fluor is in fluors of Userm
  if(!(fluor %in% fluors)){
    stop(paste0("Error: the input fluor: ",fluor," is not in the Userm."))
  }
  #check scale
  if(!(scale %in% c("Linear","Log10","Arcsinh"))){
    stop(paste0("Error: the input scale should be one of 'Linear','Log10','Arcsinh'."))
  }
  #check min and max if is numeric
  if (!is.numeric(min) || !is.numeric(max)) {
    stop("Error: 'min' and 'max' must be numeric.")
  }
  #check if min is smaller than max
  if (min >= max) {
    stop("Error: 'min' must be smaller than 'max'.")
  }
  #check if cofactor is numeric and positive
  if (!is.numeric(cofactor) || cofactor <= 0) {
    stop("Error: 'cofactor' must be a positive numeric value.")
  }

  Userm$Scale_df[fluor,"min"] = min
  Userm$Scale_df[fluor,"max"] = max
  Userm$Scale_df[fluor,"scale"] = scale
  Userm$Scale_df[fluor,"cofactor"] = cofactor

  return(Userm)

}
