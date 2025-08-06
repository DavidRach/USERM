#' Query and summarize Signature list information
#'
#' This function reads a serialized RDS file containing a list of Signature objects,
#' extracts key metadata from each Signature, and returns a structured data frame.
#'
#' @return A `data.frame` with the following columns:
#' \describe{
#'   \item{id}{Unique identifier of the Signature}
#'   \item{PrimaryName}{Primary name of the Signature}
#'   \item{SecondaryName}{Secondary name of the Signature}
#'   \item{detectors}{Number of detectors in the Signature}
#'   \item{instrument}{Instrument associated with the Signature}
#' }
#'
#' @examples
#' \dontrun{
#'   sig_info <- querySig()
#'   head(sig_info)
#' }
#'
#' @export

querySig = function(){

  Sig_list = readRDS(system.file("sig", "Sig_list.rds", package = "USERM"))

  info_df = as.data.frame(matrix(nrow = length(Sig_list), ncol = 5))
  colnames(info_df) = c("id","PrimaryName","SecondaryName","detectors","instrument")

  for (i in 1:length(Sig_list)) {
    info_df[i,"id"] = Sig_list[[i]]$id
    info_df[i,"PrimaryName"] = Sig_list[[i]]$PrimaryName
    info_df[i,"SecondaryName"] = Sig_list[[i]]$SecondaryName
    info_df[i,"detectors"] = length(Sig_list[[i]]$Signature)
    info_df[i,"instrument"] = Sig_list[[i]]$instrument
  }

  return(info_df)
}
