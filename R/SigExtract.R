
#' Extract a normalized signature from a single-stained control (SCC)
#'
#' This function computes a normalized signature vector by comparing positive and negative control data
#' for a set of specified detector channels. It applies a summary method (e.g., \code{median}, \code{mean})
#' to each column in the positive and negative data frames, calculates the difference, and normalizes the result.
#'
#' @param df_pos A data frame containing positive control signal data (cells x detectors).
#' @param df_neg A data frame containing negative control signal data (cells x detectors).
#' @param cols A character vector specifying the column names (detector channels) to be used for signature extraction.
#' @param method A function used to summarize each column (e.g., \code{median}, \code{mean}).
#' @param PrimaryName A character string specifying the primary name of the fluorophore.
#' @param SecondaryName A character string specifying the secondary name of the fluorophore or marker.
#' @param id A character string used to identify the signature object.
#' @param instrument A character string used to specify the instrument used to generate the raw data, e.g. Aurora5L, Xenith.
#'
#' @return A list containing:
#' \item{id}{The identifier string.}
#' \item{PrimaryName}{The primary name of the fluorophore.}
#' \item{SecondaryName}{The secondary name of the fluorophore or marker.}
#' \item{Signature}{A normalized numeric vector representing the extracted signature.}
#'
#' @examples
#' \dontrun{
#' df_pos <- data.frame(D1 = rnorm(100, 100), D2 = rnorm(100, 200))
#' df_neg <- data.frame(D1 = rnorm(100, 10), D2 = rnorm(100, 20))
#' sig <- ExtractSig(df_pos, df_neg, cols = c("D1", "D2"), method = mean,
#'                   PrimaryName = "FITC", SecondaryName = "CD3", id = "FITC_CD3")
#' }
#' @export

ExtractSig <- function(df_pos, df_neg, cols, method, PrimaryName, SecondaryName, id, instrument){

  missing_cols <- setdiff(cols, colnames(df_pos))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following required columns are missing from the df_pos:",
               paste(missing_cols, collapse = ", ")))
  }
  missing_cols <- setdiff(cols, colnames(df_neg))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following required columns are missing from the df_pos:",
               paste(missing_cols, collapse = ", ")))
  }

  sig_pos = apply(df_pos[,cols], 2, method, na.rm = TRUE)
  sig_neg = apply(df_neg[,cols], 2, method, na.rm = TRUE)

  sig_dif = sig_pos - sig_neg

  sig = sig_dif/max(sig_dif)

  return(list(id = id,
              PrimaryName = PrimaryName,
              SecondaryName = SecondaryName,
              Signature = sig,
              instrument = instrument))
}
