
#' Estimate Hotspot Matrix
#'
#' This function estimates a hotspot matrix based on the user's input matrix `A`, using selected detectors and fluorophores.
#' It computes the Moore-Penrose pseudoinverse of the submatrix defined by `detectors` and `fluors`, then applies a square root transformation to the absolute values.
#'
#' @param Userm A list containing the following elements:
#'   \itemize{
#'     \item \code{detectors}: A character vector specifying the row names (detectors) to extract from matrix \code{A}.
#'     \item \code{fluors}: A character vector specifying the column names (fluorophores) to extract from matrix \code{A}.
#'     \item \code{A}: A numeric matrix with row names and column names corresponding to detectors and fluorophores.
#'   }
#'
#' @return A numeric matrix representing the estimated hotspot matrix. The matrix is derived from the square root of the absolute values of the pseudoinverse of the selected submatrix.
#'
#' @importFrom MASS ginv
#' @export

EstimateHotspotMtx = function(Userm){

  detectors = Userm$detectors
  fluors = Userm$fluors
  A = Userm$A
  A = A[detectors, fluors, drop = FALSE]
  A_pinv = ginv(A)
  colnames(A_pinv) = rownames(A)
  rownames(A_pinv) = colnames(A)
  H_mtx = sqrt(abs(A_pinv))

  return(H_mtx)
}
