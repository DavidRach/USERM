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
#' @return A list representing the initialized residual object, containing the input matrices and placeholders for analysis results.
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

CreateRes = function(id, R, A){

  #check if id is provided
  if (missing(id)) {
    message("Note: Parameter 'id' was not provided. Using default value 'x'.")
    id = "x"
  }

  #check if A has good colnames and rownames
  detectors = rownames(A)
  fluors = colnames(A)
  if(is.null(detectors)){
    stop("Error: A (detectors x fluors) has no rownames, cannot check if A is compatible with R.")
  }
  if(is.null(fluors)){
    message("Warning: A (detectors x fluors) has no colnames, defalt fluor names are assigned.")
    colnames(A) = paste0("f",c(1:ncol(A)))
  }

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
  R = R[,rownames(A)]

  #create Res object
  Res = list(id = id,
                R = R,
                A = A,
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
