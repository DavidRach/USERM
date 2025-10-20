#' Estimate Hotspot Matrix
#'
#' Computes a hotspot matrix based on a user-defined matrix `A`, using specified detectors and fluorophores.
#' The function extracts a submatrix from `A`, calculates its Gram matrix, applies the Moore-Penrose pseudoinverse,
#' and returns the square root of the absolute values of the result.
#'
#' @param A A numeric matrix with named rows and columns representing detectors and fluorophores.
#'
#' @return A numeric matrix representing the estimated hotspot matrix, derived from the pseudoinverse of the Gram matrix of the selected submatrix.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Computes the Gram matrix: \code{t(A) \%*\% A}.
#'   \item Applies the Moore-Penrose pseudoinverse to the Gram matrix.
#'   \item Takes the square root of the absolute values of the pseudoinverse matrix.
#' }
#'
#' @importFrom MASS ginv
#' @export

EstimateHotspotMtx = function(A){

  A_pinv = ginv(A)
  colnames(A_pinv) = rownames(A)
  rownames(A_pinv) = colnames(A)
  H_mtx = (t(A) %*% A)
  H_mtx = ginv(H_mtx)
  H_mtx = sqrt(abs(H_mtx))
  colnames(H_mtx) = colnames(A)
  rownames(H_mtx) = colnames(A)
  return(H_mtx)
}
