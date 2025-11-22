#' Rename a fluorophore in a Userm object
#'
#' This function renames a fluorophore across all relevant components of a
#' `Userm` object, including `fluors`, matrices, data frames, and result lists.
#' It ensures consistency of naming throughout the object and updates the
#' `Rename_table` accordingly.
#'
#' @param Userm A `Userm` object containing fluorophore information, matrices,
#'   and results. Must include elements `fluors`, `A`, `Scale_df`,
#'   `Intensity_mtx`, `Res`, and `Rename_table`.
#' @param raw_name A character string specifying the existing fluorophore name
#'   to be replaced. Must exist in `Userm$fluors` or `Userm$Rename_table`.
#' @param new_name A character string specifying the new fluorophore name.
#'   Must not already exist in `Userm$fluors`.
#'
#' @return The updated `Userm` object with the fluorophore renamed consistently
#'   across all relevant slots.
#'
#' @details
#' The function performs several checks:
#' - Ensures all arguments are provided.
#' - Validates that `raw_name` exists in `Userm$fluors`.
#' - Ensures `new_name` does not already exist in `Userm$fluors`.
#' - Updates names in `fluors`, `A`, `Scale_df`, `Intensity_mtx`, and `Res`.
#' - Updates the `Rename_table` entry for the fluorophore.
#'
#' If the renaming is successful, a message is printed to the console.
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' }
#'
#' @export

RenameFluor <- function(Userm, raw_name, new_name) {
  options(error = NULL)

  # Check if all arguments are provided
  if (missing(Userm) || missing(raw_name) || missing(new_name)) {
    stop("You must provide Userm, raw_name, and new_name")
  }

  # Check if raw_name exists in Userm$fluors
  if (!(raw_name %in% Userm$fluors)) {
    stop(paste("The raw_name", raw_name, "is not found in Userm$fluors"))
  }

  # Check if new_name already exists in Userm$fluors
  if (new_name %in% Userm$fluors) {
    stop(paste("The new_name", new_name, "already exists in Userm$fluors"))
  }

  # Perform the renaming
  Userm$fluors[Userm$fluors == raw_name] = new_name
  colnames(Userm$A)[colnames(Userm$A) == raw_name] = new_name
  rownames(Userm$Scale_df)[rownames(Userm$Scale_df) == raw_name] = new_name
  rownames(Userm$Intensity_mtx)[rownames(Userm$Intensity_mtx) == raw_name] = new_name

  Userm$Res[[raw_name]]$id = new_name
  Userm$Res[[raw_name]]$fluors = new_name
  colnames(Userm$Res[[raw_name]]$A) = new_name
  names(Userm$Res)[names(Userm$Res) == raw_name] = new_name

  rename_indice = grep(raw_name,Userm[["Rename_table"]]$raw_name)
  if(length(rename_indice) == 0){
    rename_indice = grep(raw_name,Userm[["Rename_table"]]$new_name)
  }
  if(length(rename_indice) == 0){
    stop(paste("The raw_name", raw_name, "is not found in Userm$Rename_table"))
  }else{
    Userm[["Rename_table"]]$new_name[rename_indice] = new_name
  }

  # Optionally, print a message
  message(paste("Successfully renamed", raw_name, "as", new_name))

  # Return the updated object
  return(Userm)
}
