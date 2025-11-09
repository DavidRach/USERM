
#' Estimate Spillover Spreading Matrix (SSM) from Unmixing and Control Data
#'
#' Computes a full Spillover Spreading Matrix (SSM) for a set of fluorophores using precomputed
#' positive/negative control data and an unmixing matrix. This matrix quantifies the extent of signal
#' spreading from each fluorophore into others, aiding in panel design and instrument evaluation
#' in spectral flow cytometry.
#'
#' @param SSM_fluor A character vector of fluorophore names to include in the SSM. Each must have a corresponding `SSMObj_*.rds` file.
#' @param A A numeric unmixing matrix (detectors × fluorophores) with row names as detector channels and column names as fluorophores.
#' @param custom_ssm_dir Optional path to a directory containing custom `SSMObj_*.rds` files. If `NULL`, only built-in files from the `USERM` package are used.
#'
#' @return A numeric matrix of spillover spreading scores, with rows and columns named by `SSM_fluor`.
#' Each entry `ssm[i, j]` represents the spreading from fluorophore `i` into fluorophore `j`.
#' Diagonal entries are set to 0.
#'
#' @details
#' For each fluorophore in `SSM_fluor`, the function loads its corresponding SSM object containing
#' positive and negative control data. It performs detector-to-fluorophore unmixing using the pseudoinverse
#' of matrix `A`, then calculates quantile-based spread metrics for each fluorophore pair.
#'
#' The spillover spreading score is computed as:
#' \deqn{SS = \sqrt{(S_{\sigma}^2 - R_{\sigma}^2) / \Delta F}}
#' where \eqn{S_{\sigma}} and \eqn{R_{\sigma}} are the 84th–50th percentile spreads in positive and negative samples,
#' and \eqn{\Delta F} is the median signal increase in the primary fluorophore.
#'
#' If \eqn{\Delta F = 0}, the score is set to `NA`. Negative intermediate values yield imaginary scores,
#' which are returned as negative square roots.
#'
#' @importFrom MASS ginv
#' @examples
#' \dontrun{
#' A <- matrix(runif(25), nrow = 5)
#' rownames(A) <- paste0("Det", 1:5)
#' colnames(A) <- c("FITC", "PE", "APC", "PerCP", "BV421")
#' EstimateSSM(SSM_fluor = c("FITC", "PE", "APC"), A = A, custom_ssm_dir = "~/custom_ssm")
#' }
#' @export


EstimateSSM = function(SSM_fluor,A,custom_ssm_dir = NULL){
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



  #create empty ssm
  ssm = matrix(nrow = length(SSM_fluor), ncol = length(SSM_fluor), dimnames = list(SSM_fluor,SSM_fluor))


  for (fluor in SSM_fluor) {
    # fluor = SSM_fluor[2]
    print(fluor)
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

    #calculate ss
    for(fluor_negchannel in SSM_fluor){
      # fluor_negchannel = SSM_fluor[1]
      if(fluor_negchannel == fluor){
        ssm[fluor,fluor_negchannel] = 0
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
          ss = NA
        }else{
          ss_2 = delta_sigma_F_neg_2 / delta_F_pos
          if(ss_2 < 0){
            ss = - sqrt(-ss_2)
          }else{
            ss = sqrt(ss_2)
          }
        }

        ssm[fluor,fluor_negchannel] = ss

      }

    }

  }

  return(ssm)
}

