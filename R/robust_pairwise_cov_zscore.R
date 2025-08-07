#' Robust Pairwise Covariance Matrix Using Z-score Outlier Removal
#'
#' Computes a robust covariance matrix by calculating pairwise covariances between columns,
#' excluding outliers identified via Z-score thresholding. For each pair of variables,
#' rows containing outliers in either variable are excluded from the covariance calculation.
#' Diagonal entries are filled with the variance of each column after excluding its own outliers.
#'
#' @param data A numeric matrix or data frame where rows are observations and columns are variables.
#' @param z_thresh A numeric value specifying the Z-score threshold for outlier detection.
#'        Observations with absolute Z-scores greater than this value are considered outliers.
#'        Default is 3.
#'
#' @return A symmetric covariance matrix with robust pairwise estimates.
#'
#' @examples
#' set.seed(123)
#' data <- matrix(rnorm(100 * 5), ncol = 5)
#' data[1:5, ] <- data[1:5, ] + rnorm(5 * 5, sd = 10)  # Inject outliers
#' cov_matrix <- robust_pairwise_cov_zscore(data)
#' print(cov_matrix)
#'
#' @export

robust_pairwise_cov_zscore <- function(data, z_thresh = 3) {
  data <- as.data.frame(data)
  n <- ncol(data)
  cov_matrix <- matrix(NA, nrow = n, ncol = n)
  colnames(cov_matrix) <- colnames(data)
  rownames(cov_matrix) <- colnames(data)

  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      x <- data[[i]]
      y <- data[[j]]

      # cal Z-score
      z_x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
      z_y <- (y - mean(y, na.rm = TRUE)) / sd(y, na.rm = TRUE)

      # label outlier
      x_outlier <- abs(z_x) > z_thresh
      y_outlier <- abs(z_y) > z_thresh

      # exclude outlier
      keep <- !(x_outlier | y_outlier)

      # cal cov
      cov_ij <- cov(x[keep], y[keep], use = "complete.obs")
      cov_matrix[i, j] <- cov_ij
      cov_matrix[j, i] <- cov_ij  # Symmetric filling
    }
  }

  # Fill the diagonal
  for (i in 1:n) {
    x <- data[[i]]
    z_x <- (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
    x_outlier <- abs(z_x) > z_thresh
    keep <- !x_outlier
    cov_matrix[i, i] <- var(x[keep], na.rm = TRUE)
  }

  return(cov_matrix)
}
