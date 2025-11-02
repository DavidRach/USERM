#' Estimate slope and intercept matrices from residuals
#'
#' This function performs robust estimation of slope and intercept matrices based on residuals from a signal reconstruction model.
#' It bins the data using either absolute or percentile-based methods, computes robust covariance matrices for each bin,
#' and fits robust linear models to estimate the relationship between covariance values and bin midpoints.
#'
#' @param Res A result object created by \code{\link{CreateRes}}, containing raw signal matrix \code{R}, unmixing matrix \code{A}, and detector information.
#' @param count_thre An integer specifying the minimum number of data points required in a bin to perform covariance estimation.
#' @param bin_num An integer specifying the number of bins to divide the data into. Default is 30.
#' @param bin_method A character string specifying the binning method. Must be either \code{"absolute"} or \code{"percentile"}. Default is \code{"percentile"}.
#' @param z_thre A numeric threshold for removing outliers based on z-scores of the first column of matrix \code{B}. Default is 3.
#' @param ... Additional arguments passed to internal functions (currently unused).
#'
#' @return The updated \code{Res} object with the following fields added:
#' \item{bin_mids}{Midpoints of each bin used for slope estimation.}
#' \item{bin_counts}{Number of data points in each bin.}
#' \item{cov_matrices}{Robust covariance matrices computed for each bin.}
#' \item{interceptMtx}{Matrix of intercepts estimated from robust linear models.}
#' \item{slopMtx}{Matrix of slopes estimated from robust linear models.}
#' \item{par}{A list containing binning parameters used in the estimation.}
#'
#' @examples
#' \dontrun{
#' Res <- CreateRes("F1", R, A)
#' Res <- SlopEstimation(Res, count_thre = 30, bin_num = 30, bin_method = "percentile")
#' }
#' @export
#' @importFrom MASS ginv rlm cov.rob
#' @importFrom stats quantile cov
#' @importFrom progress progress_bar

SlopEstimation = function(Res, count_thre, bin_num = 30, bin_method = "percentile", z_thre = 3, ...){
  set.seed(123)

  #check if R and A are available
  #may skip this

  #calculate Residual and B
  R = Res$R
  A = Res$A #A (detectors x fluors)
  A_AF = Res$A_AF #A (detectors x fluors)
  A = cbind(A,A_AF)
  detectors = Res$detectors

  #unmix
  R = t(R) #R (detectors x cells)
  A_pinv = ginv(A)
  B =  A_pinv %*% R #B (fluors x cells)

  #for following analysis, transform B and R
  R = t(R)
  B = t(B)
  R_explained = B %*% t(A)
  Residual = R - R_explained

  B = B[, 1, drop = FALSE]

  #set default for count_thre
  if (missing(count_thre)) {
    count_thre = dim(A)[1]
  }

  #remove outliers based on B
  z_scores = scale(B[,1])
  mask = abs(z_scores) <= z_thre
  # table(mask)
  B = B[mask, , drop = FALSE]
  Residual = Residual[mask, , drop = FALSE]

  #calculate bin_size
  if(bin_method == "absolute"){
    bin_size = (max(B) - min(B)) / bin_num
  }else if(bin_method == "percentile"){
    bin_size = 1 / bin_num
  }else{
    stop("Error: bin_method should be 'absolute' or 'percentile'.")
  }

  #initial cov_matrices and bin_mids
  message("calculating cov matrix...")
  cov_matrices = array(0, dim = c(dim(A)[1], dim(A)[1], bin_num)) #(51,51,bin_num)
  bin_mids = c()
  bin_counts = c()
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                         total = bin_num,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 100)
  #loop bin and fill cov_matrices and bin_mids
  for (i in 1:bin_num) {
    pb$tick()
    # i = 2
    if(bin_method == "absolute"){
      bin_min = min(B) + (i - 1) * bin_size  - 1 * bin_size
      bin_max = bin_min + bin_size + 1 * bin_size
      bin_mid = (bin_min + bin_max) / 2

    }else if(bin_method == "percentile"){
      bin_min = quantile(B,(i - 1) * bin_size)[[1]]
      bin_max = quantile(B,(i) * bin_size)[[1]]
      bin_mid = quantile(B,(i - 0.5) * bin_size)[[1]]
    }

    bin_mask = (B >= bin_min) & (B < bin_max)
    bin_points = B[bin_mask]

    bin_counts = c(bin_counts, length(bin_points))

    if (length(bin_points) > count_thre) {
      residuals_in_bin = Residual[bin_mask,]
      dim(residuals_in_bin)
      # cov_matrices[ , ,i] = MASS::cov.rob(residuals_in_bin)$cov #will cause sigularity issue when cal inv matrix of cov matrix
      cov_matrices[ , ,i] = robust_pairwise_cov_zscore(residuals_in_bin,z_thresh = z_thre)
      bin_mids = c(bin_mids, bin_mid)
    }else{
      bin_mids = c(bin_mids, NA)
    }
  }
  # save bin_mids, bin_counts, cov_matrices,
  Res$bin_mids = bin_mids
  Res$bin_counts = bin_counts
  Res$cov_matrices = cov_matrices

  if(any(is.na(bin_mids))){
    message(paste0("There are ",sum(is.na(bin_mids)), " bins do not have enough sample for estiamtion."))
  }
  cov_matrices = cov_matrices[,,(!is.na(bin_mids))]
  bin_mids = bin_mids[!is.na(bin_mids)]

  #initial interceptMtx and slopMtx
  message("calculating slop matrix...")
  interceptMtx = array(0, dim = c(dim(A)[1], dim(A)[1]),dimnames = list(detectors,detectors))#(51,51)
  slopMtx = interceptMtx #(51,51)

  # loop detector pairs to fill interceptMtx and slopMtx
  for (i in 1:dim(A)[1]) {
    for (j in 1:dim(A)[1]) {

      covs = cov_matrices[i,j,]
      model_robust = rlm(covs ~ bin_mids)
      interceptMtx[i,j] = coef(model_robust)[1]
      slopMtx[i,j] = coef(model_robust)[2]
    }
  }
  # save bin_num, count_thre, interceptMtx and slopMtx into Res and return Res
  Res$par$bin_num = bin_num
  Res$par$bin_method = bin_method
  Res$par$count_thre = count_thre
  Res$interceptMtx = interceptMtx
  Res$slopMtx = slopMtx

  return(Res)
}

