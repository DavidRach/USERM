
#' Estimate Cosine Similarity Matrix for Fluorochromes
#'
#' This function computes the pairwise cosine similarity between fluorochrome spectral profiles
#' based on a user-defined matrix of detector responses. It extracts a submatrix from the full
#' spectral matrix using specified detectors and fluorochromes, and returns a symmetric similarity
#' matrix where each entry represents the cosine similarity between two fluorochromes.
#'
#' @param A A numeric matrix with named rows and columns representing detectors and fluorophores.
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
EstimateSimilarityMtx = function(A){

  cos_sim_matrix <- cosine(A)
  colnames(cos_sim_matrix) = colnames(A)
  rownames(cos_sim_matrix) = colnames(A)
  return(cos_sim_matrix)
}
