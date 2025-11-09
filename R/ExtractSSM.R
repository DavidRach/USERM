
#' Extract Spillover Spreading Matrix (SSM) related object from Positive and Negative Controls
#'
#' Computes a Spillover Spreading Matrix (SSM) related object by applying a user-defined
#' summary method to selected signal channels from positive and negative control datasets.
#'
#' @param df_pos A data frame containing positive control samples (e.g., single-stained beads).
#' @param df_neg A data frame containing negative control samples (e.g., unstained or background).
#' @param cols A character vector specifying the signal channels (column names) to include in SSM calculation.
#' @param method A summary function (e.g., `mean`, `median`) applied to each channel. Must accept `na.rm = TRUE`.
#' @param PrimaryName A character string specifying the primary name of the fluorophore.
#' @param SecondaryName A character string specifying the secondary name of the fluorophore or marker.
#' @param id A unique identifier for the SSM signature.
#' @param instrument A character string describing the cytometry instrument used.
#' @param Source A character string indicating the source or origin of the dataset.
#' @param Note Optional notes or annotations. Default is `NA`.
#'
#' @return A list containing:
#' \item{id}{The provided identifier.}
#' \item{PrimaryName}{Label for the primary condition.}
#' \item{SecondaryName}{Label for the secondary condition.}
#' \item{Signature}{A numeric vector of normalized spreading scores across channels.}
#' \item{instrument}{Instrument description.}
#' \item{Source}{Data source.}
#' \item{Note}{Optional notes.}
#' \item{df_pos}{Sampled positive control data (max 1000 rows).}
#' \item{df_neg}{Sampled negative control data (max 1000 rows).}
#'
#' @details
#' The Spillover Spreading Matrix (SSM) quantifies the extent of signal spreading from one channel into others,
#' which is critical for panel design and instrument performance evaluation in spectral flow cytometry.
#'
#' Input data frames are subsetted to the specified channels and randomly sampled (up to 1000 rows)
#' to reduce output size and facilitate downstream visualization.
#'
#' @importFrom dplyr sample_n
#' @examples
#' \dontrun{
#' pos <- data.frame(CD3 = rnorm(2000), CD19 = rnorm(2000))
#' neg <- data.frame(CD3 = rnorm(2000), CD19 = rnorm(2000))
#' ExtractSSM(df_pos = pos, df_neg = neg, cols = c("CD3", "CD19"),
#'            method = median, PrimaryName = "FITC", SecondaryName = "CD3",
#'            id = "FITC_CD3", instrument = "Aurora", Source = "PanelX")
#' }
#' @export
#' @importFrom dplyr sample_n


ExtractSSM <- function(df_pos, df_neg, cols, method, PrimaryName, SecondaryName, id,
                       instrument, Source, Note = NA){

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

  #sample df
  df_pos = df_pos[,cols]
  set.seed(111)
  if (nrow(df_pos)>1000) {
    df_pos = dplyr::sample_n(df_pos,size = 1000,replace = F)
  }

  df_neg = df_neg[,cols]
  set.seed(111)
  if (nrow(df_neg)>1000) {
    df_neg = dplyr::sample_n(df_neg,size = 1000,replace = F)
  }

  return(list(id = id,
              PrimaryName = PrimaryName,
              SecondaryName = SecondaryName,
              Signature = sig,
              instrument = instrument,
              Source = Source,
              Note = Note,
              df_pos = df_pos,
              df_neg = df_neg))
}
