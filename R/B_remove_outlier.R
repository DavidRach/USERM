#' Remove outliers from a data frame column using a MAD-based robust filter
#'
#' This function removes outliers from a specified column of a data frame
#' using a robust threshold based on the Median Absolute Deviation (MAD).
#' Values outside the interval \code{median ± 3 * MAD} are excluded.
#' The method is resistant to extreme values and does not assume normality.
#'
#' @param B A data frame containing the data to be filtered.
#' @param col A character string specifying the column name in \code{B}
#'   on which outlier removal should be performed.
#'
#' @return A filtered data frame with outliers removed from the specified column.
#'
#' @details
#' The MAD-based rule is a robust alternative to standard deviation filtering.
#' It is particularly suitable for skewed or heavy-tailed distributions,
#' where traditional z-score filtering may fail. If the column contains
#' no variability (MAD = 0), the function will return the original data frame.
#'
#' @examples
#' df <- data.frame(value = c(1, 2, 3, 100))
#' B_remove_outlier(df, "value")
#'
#' @export

B_remove_outlier = function(B, col, k = 5){
  x <- B[,col]

  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, constant = 1, na.rm = TRUE)

  lower <- med - 5 * mad_val
  upper <- med + 5 * mad_val

  B <- B[x >= lower & x <= upper, ]

  return(B)
}

