#' Retrieve Result Object by ID
#'
#' This function loads a serialized result object from the internal package directory
#' based on a given ID. The result object is expected to be stored as an RDS file
#' named \code{ResObj_<id>.rds} within the \code{res} folder of the \code{USERM} package.
#'
#' @param id A character string representing the ID of the result object to retrieve.
#'
#' @details
#' The function constructs the file path using \code{system.file()} and checks whether
#' the corresponding RDS file exists. If the file is found, it is read and returned.
#' If not, the function throws an error indicating that the result object could not be found.
#'
#' @return An R object containing the result data associated with the specified ID.
#'
#' @examples
#' \dontrun{
#'   res <- getRes("SCC_Bead_CD4_NFR700")
#'   print(res)
#' }
#'
#' @export

getRes = function(id,custom_dir = NULL){

  if(is.null(custom_dir)){
    file_addr = system.file("res", paste0("ResObj_",id,".rds"), package = "USERM")
  }else{
    file_addr = paste0(custom_dir,"/ResObj_",id,".rds")
  }


  if(file.exists(file_addr)){
    Res = readRDS(file_addr)
    return(Res)
  }else{
    stop("The Res object for the id is not found. Please check the id and try again.")
  }
}
