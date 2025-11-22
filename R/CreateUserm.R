
#' Create a userm object
#'
#' This function initializes a userm object (`Userm`) based on an unmixing matrix `A`.
#' It validates the structure of `A`, assigns default fluor names if missing, and prepares
#' placeholders for storing results, scaling parameters, and intensity values.
#'
#' @param A A numeric unmixing matrix (detectors × fluors), typically extracted from a single-stained control (SCC) FCS file.
#'           Row names should represent detector channels, and column names should represent fluorophores.
#'
#' @return A list representing the userm object, containing:
#' \itemize{
#'   \item \code{A}: the input unmixing matrix.
#'   \item \code{detectors}: vector of detector names.
#'   \item \code{fluors}: vector of fluorophore names.
#'   \item \code{Res}: a named list to store res objects for each fluor.
#'   \item \code{Scale_df}: a data frame for storing scaling parameters (min, max, scale, cofactor).
#'   \item \code{Intensity_mtx}: a data frame for storing default intensity values.
#' }
#'
#' @examples
#' \dontrun{
#' A <- matrix(runif(20), nrow = 4, ncol = 5)
#' rownames(A) <- paste0("D", 1:4)
#' colnames(A) <- paste0("F", 1:5)
#' Userm <- CreateUserm(A)
#' }
#' @export

CreateUserm = function(A){

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

  Res_list = list()
  for (i in 1:length(fluors)) {
    Res_list[[i]] = list()
  }
  names(Res_list) = fluors

  Intensity_mtx = as.data.frame(matrix(data = 100,
                                            nrow = length(fluors),
                                            ncol = 1,
                                            dimnames = list(fluors,c("P1"))))

  Scale_df = as.data.frame(matrix(data = NA,
                                  nrow = length(fluors),
                                  ncol = 4,
                                  dimnames = list(fluors,c("min","max","scale","cofactor"))))
  Scale_df$min = 0
  Scale_df$max = 1000
  Scale_df$scale = "Linear"
  cofactor = 10

  Rename_table = as.data.frame(matrix(data = NA,
                                      nrow = length(fluors),
                                      ncol = 2,
                                      dimnames = list(fluors,c("raw_name","new_name"))))
  Rename_table$raw_name = rownames(Rename_table)
  Rename_table$new_name = rownames(Rename_table)

  Userm = list(A = A,
                 detectors = detectors,
                 fluors = fluors,
                 Res = Res_list,
                 Scale_df = Scale_df,
                 Intensity_mtx = Intensity_mtx,
                 Rename_table = Rename_table)

  return(Userm)
}
