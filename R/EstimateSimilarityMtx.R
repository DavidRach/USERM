
#' Estimate Cosine Similarity Matrix for Fluorochromes
#'
#' This function computes the pairwise cosine similarity between fluorochrome spectral profiles
#' based on a user-defined matrix of detector responses. It extracts a submatrix from the full
#' spectral matrix using specified detectors and fluorochromes, and returns a symmetric similarity
#' matrix where each entry represents the cosine similarity between two fluorochromes.
#'
#' @param Userm A list containing the following components:
#'   \itemize{
#'     \item \code{detectors}: A character vector of detector names to be used as rows.
#'     \item \code{fluors}: A character vector of fluorochrome names to be used as columns.
#'     \item \code{A}: A numeric matrix with rows representing detectors and columns representing fluorochromes.
#'   }
#'
#' @return A symmetric numeric matrix of cosine similarity values between fluorochromes.
#'   The matrix has fluorochrome names as both row and column names.
#'
#' @examples
#' \dontrun{
#' EstimateSimilarityMtx(Userm)
#' }
#'
#' @importFrom lsa cosine
#' @export
EstimateSimilarityMtx = function(Userm){
  detectors = Userm$detectors
  fluors = Userm$fluors
  A = Userm$A
  A = A[detectors, fluors, drop = FALSE]
  cos_sim_matrix <- cosine(A)
  colnames(cos_sim_matrix) = colnames(A)
  rownames(cos_sim_matrix) = colnames(A)
  return(cos_sim_matrix)
}
