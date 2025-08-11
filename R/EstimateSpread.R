
#' Estimate Signal Spread for a Given Population
#'
#' This function estimates the signal spread (variance components) for a specified population
#' using the USERM result object. It computes weighted intercept and slope variances across
#' detectors and fluorochromes based on the system matrix and signal intensities.
#'
#' @param Userm A list containing USERM analysis results. Must include:
#' \itemize{
#'   \item \code{Intensity_mtx}: A matrix of signal intensities (fluorochrome x population).
#'   \item \code{detectors}: A character vector of detector names.
#'   \item \code{fluors}: A character vector of fluorochrome names.
#'   \item \code{A}: A system matrix mapping detectors to fluorochromes.
#'   \item \code{Res}: A named list of result objects for each fluorochrome, each containing \code{interceptMtx} and \code{slopMtx}.
#' }
#' @param population_id A character string specifying the population ID to estimate spread for. Must match a column name in \code{Userm$Intensity_mtx}.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{population_id}: The input population ID.
#'   \item \code{A}: The system matrix used.
#'   \item \code{fluors}: Fluorochrome names.
#'   \item \code{detectors}: Detector names.
#'   \item \code{intensity_matrix}: Signal intensity values for the specified population.
#'   \item \code{intercept_sigma2_col}: Estimated intercept variance per channel.
#'   \item \code{slop_sigma2_col}: Estimated slope variance per channel.
#' }
#'
#' @examples
#' \dontrun{
#'   EstimateSpread(Userm, "Population_1")
#' }
#'
#' @importFrom MASS ginv
#' @export

EstimateSpread = function(Userm,population_id){
  if (!population_id %in% colnames(Userm$Intensity_mtx)) {
    stop(paste0("The ",population_id," column is missing in Userm$Intensity_mtx."))
  }

  detectors = Userm$detectors
  fluors = Userm$fluors
  A = Userm$A
  A = A[detectors, fluors, drop = FALSE]
  intensity_matrix = Userm$Intensity_mtx[,population_id,drop = FALSE]

  A_pinv = ginv(A)
  colnames(A_pinv) = rownames(A)
  rownames(A_pinv) = colnames(A)

  #calculate weighted intercept
  weighted_matrix <- array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)
  for (i in 1:ncol(A)) {
    fluor = colnames(A)[i]
    intercept_matrix = Userm$Res[[fluor]]$interceptMtx
    intercept_matrix = intercept_matrix[detectors,detectors]
    weighted_matrix[,i] = diag((A_pinv %*% intercept_matrix) %*% t(A_pinv))
  }
  intercept_sigma2_col = matrix(apply(weighted_matrix,1,median),ncol = 1)
  # intercept_sigma_col = sqrt(abs(intercept_sigma2_col))
  rownames(intercept_sigma2_col) = colnames(A)
  # rownames(intercept_sigma_col) = colnames(A)

  #calculate weighted slop
  weighted_matrix = array(0, dim = c(ncol(A), ncol(A))) #(channel, scc)
  for (i in 1:ncol(A)) {
    fluor = colnames(A)[i]
    slop_matrix = Userm$Res[[fluor]]$slopMtx
    slop_matrix = slop_matrix[detectors,detectors]
    weighted_matrix[,i] = diag((A_pinv %*% slop_matrix) %*% t(A_pinv)) * intensity_matrix[fluor,1]
  }
  slop_sigma2_col = matrix(apply(weighted_matrix,1,sum),ncol = 1)
  # slop_sigma_col = sqrt(abs(slop_sigma2_col))
  rownames(slop_sigma2_col) = colnames(A)
  # rownames(slop_sigma_col) = colnames(A)

  Spr = list(population_id = population_id,
             A = A,
             fluors = fluors,
             detectors = detectors,
             intensity_matrix = intensity_matrix,
             intercept_sigma2 = intercept_sigma2_col,
             slop_sigma2 = slop_sigma2_col
             )
  return(Spr)

}

