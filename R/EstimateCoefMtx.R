#' Estimate Coefficient Matrix of Slope Components in USERM
#'
#' This function computes a tangent matrix representing the diagonalized slope components
#' of signal spread across fluorochromes, based on the system matrix and slope matrices
#' from the USERM result object.
#'
#' @param Userm A list containing USERM analysis results. Must include:
#' \itemize{
#'   \item \code{detectors}: A character vector of detector names.
#'   \item \code{fluors}: A character vector of fluorochrome names.
#'   \item \code{A}: A system matrix mapping detectors to fluorochromes.
#'   \item \code{Res}: A named list of result objects for each fluorochrome, each containing \code{slopMtx}.
#' }
#'
#' @return A square numeric matrix with fluorochrome names as both row and column names.
#' Each entry represents the diagonalized slope component projected through the pseudo-inverse
#' of the system matrix.
#'
#' @details The function uses the pseudo-inverse of the system matrix \code{A} to transform
#' slope matrices from detector space to fluorochrome space. The diagonal elements of the
#' transformed matrices are extracted and assembled into the final tangent matrix.
#'
#' @examples
#' \dontrun{
#'   tangent_mtx <- EstimateCoefMtx(Userm)
#'   print(tangent_mtx)
#' }
#'
#' @importFrom MASS ginv
#' @export
EstimateCoefMtx = function(Userm,A=NULL){

  detectors = Userm$detectors
  fluors = Userm$fluors
  if(is.null(A)){
    A = Userm$A
  }
  A = A[detectors, fluors, drop = FALSE]
  A_pinv = ginv(A)
  colnames(A_pinv) = rownames(A)
  rownames(A_pinv) = colnames(A)

  #calculate weighted slop
  Coef_matrix = array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)
  for (i in 1:ncol(A)) {
    fluor = colnames(A)[i]
    slop_matrix = Userm$Res[[fluor]]$slopMtx
    slop_matrix = slop_matrix[detectors,detectors]
    Coef_matrix[,i] = diag((A_pinv %*% slop_matrix) %*% t(A_pinv))
  }
  colnames(Coef_matrix) = colnames(A)
  rownames(Coef_matrix) = colnames(A)
  Coef_matrix = t(Coef_matrix)
  # Coef_matrix = sqrt(abs(Coef_matrix))# Variance To Standard deviation

  for (row in 1:nrow(Coef_matrix)) {
    for (col in 1:ncol(Coef_matrix)) {
      if(Coef_matrix[row,col] < 0){
        Coef_matrix[row,col] = -sqrt(abs(Coef_matrix[row,col]))
      }else if(Coef_matrix[row,col] > 0){
        Coef_matrix[row,col] = sqrt(Coef_matrix[row,col])
      }
    }
  }
  return(Coef_matrix)
}






