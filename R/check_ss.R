
#' Compute Spillover Spreading Score (SS) Between Fluorophores
#'
#' This function calculates the spillover spreading score (SS) between a primary fluorophore (`f_pos`)
#' and a secondary fluorophore (`f_neg`) using precomputed SSM objects and an unmixing matrix `A`.
#' It supports both built-in and custom SSM sources, and performs detector-to-fluorophore unmixing
#' followed by quantile-based spread estimation.
#'
#' @param f_pos A character string specifying the positive fluorophore (positive channel).
#' @param f_neg A character string specifying the negative fluorophore (negative channel).
#' @param SSM_fluor A character vector of fluorophores required for SSM computation. Used for availability checks.
#' @param A A numeric matrix representing the detector-to-fluorophore unmixing matrix. Row names must be detector names; column names must be fluorophore names.
#' @param custom_ssm_dir Optional path to a directory containing custom SSM `.rds` files. Default is `NULL`.
#'
#' @return A list containing:
#' \item{f_pos}{positive fluorophore name.}
#' \item{f_neg}{negative fluorophore name.}
#' \item{B_pos}{Unmixed positive sample matrix (cells × fluorophores).}
#' \item{B_neg}{Unmixed negative sample matrix (cells × fluorophores).}
#' \item{R_F_neg_84}{84th percentile of `f_neg` in negative samples.}
#' \item{R_F_neg_50}{Median of `f_neg` in negative samples.}
#' \item{R_sigma_F_neg}{Spread of `f_neg` in negative samples.}
#' \item{S_F_neg_84}{84th percentile of `f_neg` in positive samples.}
#' \item{S_F_neg_50}{Median of `f_neg` in positive samples.}
#' \item{S_sigma_F_neg}{Spread of `f_neg` in positive samples.}
#' \item{S_F_pos_50}{Median of `f_pos` in positive samples.}
#' \item{R_F_pos_50}{Median of `f_pos` in negative samples.}
#' \item{delta_F_pos}{Difference in median signal of `f_pos` between positive and negative samples.}
#' \item{delta_sigma_F_neg_2}{Difference in squared spread of `f_neg` between positive and negative samples.}
#' \item{ss_2}{Intermediate SS score before square root transformation.}
#' \item{ss}{Final spillover spreading score.}
#'
#' @details
#' For the detailed explanation of returned items, please refer to Fig1A of Nguyen, Perfetto et al. Quantifying Spillover Spreading for Comparing Instrument Performance and Aiding in Multicolor Panel Design”. Cytometry A. 2013.
#'
#' The function loads SSM objects for the specified fluorophores, checks detector compatibility with the unmixing matrix,
#' performs matrix inversion and unmixing, and computes quantile-based spread metrics. The SS score is derived from
#' the change in spread of the negative channel (`f_neg`) relative to the signal increase in the positive channel (`f_pos`).
#'
#' If `f_pos == f_neg`, the score is set to zero and all intermediate metrics are returned as `NA`.
#'
#' @importFrom MASS ginv
#' @examples
#' \dontrun{
#' A <- matrix(runif(25), nrow = 5)
#' rownames(A) <- paste0("Det", 1:5)
#' colnames(A) <- c("FITC", "PE", "APC", "PerCP", "BV421")
#' check_ss(f_pos = "FITC", f_neg = "PE", SSM_fluor = c("FITC", "PE"),
#'          A = A, custom_ssm_dir = "~/custom_ssm")
#' }
#' @export


check_ss = function(f_pos, f_neg, SSM_fluor,A,custom_ssm_dir = NULL){

  #prepare A
  detector_A = rownames(A)
  fluor_A = colnames(A)
  A_pinv = ginv(A)

  #check if all SSM_fluor has corresponding data in ssm folder

  ## find available ssm_file from USERM
  all_ssm_file_builtin = list.files(system.file("ssm", package = "USERM"))
  all_ssm_file_builtin = sub("^SSMObj_", "", sub("\\.rds$", "", all_ssm_file_builtin))

  ## find available ssm_file from custom_ssm_dir
  if(is.null(custom_ssm_dir)){
    all_ssm_file_custom = c()
  }else{
    all_ssm_file_custom = dir(custom_ssm_dir)
    all_ssm_file_custom = sub("^SSMObj_", "", sub("\\.rds$", "", all_ssm_file_custom))
  }

  ## check SSM_fluor availability
  missing_cols <- setdiff(SSM_fluor, c(all_ssm_file_builtin,all_ssm_file_custom))
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following required fluors are missing from USERM and custom dir:",
               paste(missing_cols, collapse = ", ")))
  }

  #check if all SSM_fluor are in A

  missing_cols <- setdiff(SSM_fluor, fluor_A)
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following required fluors are missing from input A:",
               paste(missing_cols, collapse = ", ")))
  }


  #calculate ss
  fluor = f_pos
  # print(fluor)

  # load ssmobj
  if(fluor %in% all_ssm_file_builtin){
    tmp_ssm_obj = readRDS(system.file("ssm", paste0("SSMObj_",fluor,".rds"), package = "USERM"))
  }else{
    tmp_ssm_obj = readRDS(paste0(custom_ssm_dir,"/SSMObj_",fluor,".rds"))
  }

  #check df_pos and df_neg
  df_pos = tmp_ssm_obj$df_pos
  ##check columns of df_pos
  detector_pos = colnames(df_pos)
  missing_cols <- setdiff(detector_pos, detector_A)
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following detectors are missing from df_pos in the ssm obj of ",fluor ,":",
               paste(missing_cols, collapse = ", ")))
  }else{
    df_pos = df_pos[,detector_A]
  }

  df_neg = tmp_ssm_obj$df_neg
  ##check columns of df_neg
  detector_neg = colnames(df_neg)
  missing_cols <- setdiff(detector_neg, detector_A)
  if (length(missing_cols) > 0) {
    stop(paste("Error: The following detectors are missing from df_neg in the ssm obj of ",fluor ,":",
               paste(missing_cols, collapse = ", ")))
  }else{
    df_neg = df_neg[,detector_A]
  }

  # unmix
  R_pos = t(df_pos) #R (detectors x cells)
  B_pos =  A_pinv %*% R_pos #B (fluors x cells)
  B_pos = t(B_pos)
  colnames(B_pos) = fluor_A



  R_neg = t(df_neg) #R (detectors x cells)
  B_neg =  A_pinv %*% R_neg #B (fluors x cells)
  B_neg = t(B_neg)
  colnames(B_neg) = fluor_A

  #correct B_pos
  B_pos_correct = as.data.frame(B_pos[,c(f_pos,f_neg)])
  #correct B_pos_correct
  colnames(B_pos_correct) = c("X","Y")
  B_pos_correct = correct_to_horizontal(B_pos_correct)
  mean_neg = mean(B_neg[,c(f_neg)])
  mean_pos = mean(B_pos_correct$Y)
  B_pos[,c(f_neg)] = B_pos_correct$Y - (mean_pos - mean_neg)


  #calculate ss
  fluor_negchannel = f_neg

  if(fluor_negchannel == fluor){
    R_F_neg_84 = NA
    R_F_neg_50 = NA
    R_sigma_F_neg = NA
    S_F_neg_84 = NA
    S_F_neg_50 = NA
    S_sigma_F_neg = NA
    S_F_pos_50 = NA
    R_F_pos_50 = NA
    delta_F_pos = NA
    delta_sigma_F_neg_2 = NA
    ss_2 = NA
    ss = 0
  }else{
    R_F_neg_84 = quantile(B_neg[,fluor_negchannel],probs = 0.84)[[1]]
    R_F_neg_50 = quantile(B_neg[,fluor_negchannel],probs = 0.50)[[1]]
    R_sigma_F_neg = R_F_neg_84 - R_F_neg_50

    S_F_neg_84 = quantile(B_pos[,fluor_negchannel],probs = 0.84)[[1]]
    S_F_neg_50 = quantile(B_pos[,fluor_negchannel],probs = 0.50)[[1]]
    S_sigma_F_neg = S_F_neg_84 - S_F_neg_50

    S_F_pos_50 = quantile(B_pos[,fluor],probs = 0.50)[[1]]
    R_F_pos_50 = quantile(B_neg[,fluor],probs = 0.50)[[1]]
    delta_F_pos = S_F_pos_50 - R_F_pos_50

    delta_sigma_F_neg_2 = S_sigma_F_neg^2 - R_sigma_F_neg^2

    if(delta_F_pos == 0){
      ss_2 = NA
      ss = NA
    }else{
      ss_2 = delta_sigma_F_neg_2 / delta_F_pos
      if(ss_2 < 0){
        ss = - sqrt(-ss_2)
      }else{
        ss = sqrt(ss_2)
      }
    }
  }

  return(list(f_pos = f_pos,
              f_neg = f_neg,
              B_pos = B_pos,
              B_neg = B_neg,
              R_F_neg_84 = R_F_neg_84,
              R_F_neg_50 = R_F_neg_50,
              R_sigma_F_neg = R_sigma_F_neg,
              S_F_neg_84 = S_F_neg_84,
              S_F_neg_50 = S_F_neg_50,
              S_sigma_F_neg = S_sigma_F_neg,
              S_F_pos_50 = S_F_pos_50,
              R_F_pos_50 = R_F_pos_50,
              delta_F_pos = delta_F_pos,
              delta_sigma_F_neg_2 = delta_sigma_F_neg_2,
              ss_2 = ss_2,
              ss = ss))
}
