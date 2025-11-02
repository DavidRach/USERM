#' Create a Res object
#'
#' This function creates a structured result object (`Res`) for storing residual model information extracted from spectral flow scc file.
#' It validates the input raw signal matrix `R` and the unmixing matrix `A`, ensuring that their dimensions
#' and identifiers are compatible. If valid, it returns a list containing initialized fields for further analysis.
#'
#' @param id A character string representing the identifier for the result object. If missing, defaults to `"x"`.
#' @param R A numeric matrix of raw signal data (cells × detectors), typically extracted from an FCS file.
#' @param A A numeric unmixing matrix (detectors × fluors), typically extracted from a single-stained control (SCC) FCS file.
#'
#' @return A named list representing the initialized residual modeling object (`Res`) with the following components:
#' \describe{
#'   \item{id}{Character string identifier for the residual object. Defaults to `"x"` if not provided.}
#'   \item{R}{Numeric matrix of raw signal data (cells × detectors), subset to match `A`'s rownames.}
#'   \item{A}{Numeric unmixing matrix (detectors × fluors), with validated row and column names. For one fluorescence, here is the signature of it.}
#'   \item{detectors}{Character vector of detector names, extracted from `rownames(A)`.}
#'   \item{fluors}{Character vector of fluor names, extracted from `colnames(A)` or auto-assigned.}
#'   \item{par}{List of analysis parameters, including:
#'     \describe{
#'       \item{bin_num}{Number of bins used in residual modeling (default `NA`).}
#'       \item{bin_method}{Binning method used (e.g., `"equal-width"`, `"quantile"`; default `NA`).}
#'       \item{count_thre}{Threshold for minimum bin count (default `NA`).}
#'     }
#'   }
#'   \item{bin_mids}{Numeric vector of bin midpoints (default `NA`).}
#'   \item{bin_counts}{Numeric vector of bin counts (default `NA`).}
#'   \item{cov_matrices}{List of covariance matrices per bin (default `NA`).}
#'   \item{interceptMtx}{Matrix of intercepts from bin-wise regression (default `NA`).}
#'   \item{slopMtx}{Matrix of slopes from bin-wise regression (default `NA`).}
#' }
#' The CreateRes function is used together with SlopEstimation function to fill in most contents in the ResObj. See \code{\link{SlopEstimation}} for more detail.
#'
#' @examples
#' \dontrun{
#' R <- matrix(runif(100), nrow = 10, ncol = 5)
#' colnames(R) <- paste0("D", 1:5)
#' A <- matrix(runif(25), nrow = 5, ncol = 5)
#' rownames(A) <- paste0("D", 1:5)
#' colnames(A) <- paste0("F", 1:5)
#' Res <- CreateRes("Sample1", R, A)
#' }
#' @export

CreateRes = function(id, R, A_Target, A_AF){

  #check if id is provided
  if (missing(id)) {
    message("Note: Parameter 'id' was not provided. Using default value 'x'.")
    id = "x"
  }

  #check if A_AF has good colnames and rownames
  detectors_AF = rownames(A_AF)
  fluors_AF = colnames(A_AF)
  if(is.null(detectors_AF)){
    stop("Error: A_AF (detectors x fluors) has no rownames, cannot check if A_AF is compatible with R.")
  }
  if(is.null(fluors_AF)){
    message("Warning: A_AF (detectors x fluors) has no colnames, defalt fluor names are assigned.")
    colnames(A_AF) = paste0("f",c(1:ncol(A_AF)))
  }

  #check if A_Target has good colnames and rownames
  detectors = rownames(A_Target)
  fluors = colnames(A_Target)
  if(is.null(detectors)){
    stop("Error: A_Target (detectors x fluors) has no rownames, cannot check if A_Target is compatible with R.")
  }
  if(is.null(fluors)){
    message("Warning: A_Target (detectors x fluors) has no colnames, defalt fluor names are assigned.")
    colnames(A_Target) = paste0("f",c(1:ncol(A_Target)))
  }

  #check if A_AF has matched colnames
  missing_cols <- setdiff(detectors, detectors_AF)
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following required columns in A are missing from the A_AF:",
               paste(missing_cols, collapse = ", ")))
  }
  extra_cols <- setdiff(detectors_AF,detectors)
  if(length(extra_cols) > 0) {
    warning(paste("Warning: The following columns in A_AF are not used and will be removed:",
                  paste(extra_cols, collapse = ", ")))
  }
  A_AF = A_AF[rownames(A_Target),, drop = FALSE]

  #check if R has matched colnames
  missing_cols <- setdiff(detectors, colnames(R))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following required columns in A are missing from the R:",
               paste(missing_cols, collapse = ", ")))
  }
  extra_cols <- setdiff(colnames(R),detectors)
  if(length(extra_cols) > 0) {
    warning(paste("Warning: The following columns in R are not used and will be removed:",
                  paste(extra_cols, collapse = ", ")))
  }
  R = R[,rownames(A_Target)]

  #create Res object
  Res = list(id = id,
                R = R,
                A = A_Target,
                A_AF = A_AF,
                detectors = detectors,
                fluors = fluors,
                par = list(bin_num = NA,
                           bin_method = NA,
                           count_thre = NA),
                bin_mids = NA,
                bin_counts = NA,
                cov_matrices = NA,
                interceptMtx = NA,
                slopMtx = NA)

  return(Res)
}
