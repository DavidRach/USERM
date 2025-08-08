#' Scale Z Values with Cumulative Threshold Filtering
#'
#' This function processes a prediction grid by filtering and scaling its \code{z} values.
#' It first sorts the grid by descending \code{z}, then cumulatively sums the values until a specified threshold is reached.
#' All remaining \code{z} values beyond the threshold are set to zero.
#' Finally, the remaining \code{z} values are linearly scaled to the range [0, 1], with the maximum value becoming 1.
#' The original row order of the input grid is preserved in the output.
#'
#' @param grid A data frame containing columns \code{x}, \code{y}, and \code{z}.
#' @param threshold A numeric value between 0 and 1 specifying the cumulative proportion of \code{z} values to retain. Default is \code{0.95}.
#'
#' @return A data frame with the same structure and row order as the input \code{grid}, with filtered and scaled \code{z} values.
#'
#' @examples
#' \dontrun{
#' grid <- expand.grid(x = seq(0, 1000, length.out = 50),
#'                     y = seq(0, 1000, length.out = 50))
#' grid$z <- with(grid, exp(-((x - 500)^2 + (y - 500)^2) / 1e5))
#' scaled_grid <- scale_z(grid, threshold = 0.95)
#' }
#' @export

scale_z <- function(grid, threshold = 0.95) {
  # Save the original row order
  grid$original_order <- seq_len(nrow(grid))

  # Sort by z value in descending order
  grid_sorted <- grid[order(-grid$z), ]

  # Compute cumulative sum of z
  cum_z <- cumsum(grid_sorted$z)
  total_z <- sum(grid_sorted$z)

  # Find the index where cumulative sum exceeds the threshold
  cutoff_index <- which(cum_z / total_z > threshold)[1]

  # Set remaining z values to 0
  if (!is.na(cutoff_index)) {
    grid_sorted$z[(cutoff_index + 1):nrow(grid_sorted)] <- 0
  }

  # Linearly scale z values to [0, 1]
  max_z <- max(grid_sorted$z)
  if (max_z > 0) {
    grid_sorted$z <- grid_sorted$z / max_z
  } else {
    grid_sorted$z <- 0
  }

  # Restore original row order
  grid_restored <- grid_sorted[order(grid_sorted$original_order), ]
  grid_restored$original_order <- NULL  # Remove helper column

  return(grid_restored)
}
