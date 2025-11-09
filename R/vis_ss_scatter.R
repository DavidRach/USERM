
#' Visualize Spillover Spreading Between Fluorophores
#'
#' Generates a scatter plot comparing the signal intensities of a positive fluorophore (`f_pos`)
#' and a negative fluorophore (`f_neg`) across positive and negative control samples.
#' This visualization helps assess the extent of spillover spreading in spectral flow cytometry.
#'
#' @param ss_output A list returned by [check_ss()], containing unmixed matrices `B_pos` and `B_neg`,
#' fluorophore names `f_pos` and `f_neg`, and other spread metrics.
#'
#' @return A `ggplot` object showing a scatter plot of `f_pos` vs `f_neg` signal intensities,
#' with points colored by sample group ("Stain" for positive, "Reference" for negative).
#'
#' @details
#' The function extracts the relevant fluorophore channels from the unmixed matrices,
#' combines them into a single data frame, and plots the signal distribution using `ggplot2`.
#' This plot provides an intuitive view of how signal from `f_pos` may spread into `f_neg`,
#' aiding in panel design and spillover diagnostics.
#'
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab theme_bw
#' @examples
#' \dontrun{
#' ss <- check_ss(f_pos = "FITC", f_neg = "PE", SSM_fluor = c("FITC", "PE"),
#'                A = A, custom_ssm_dir = "~/custom_ssm")
#' vis_ss_scatter(ss)
#' }
#' @export


vis_ss_scatter = function(ss_output){

  plot_df_pos = as.data.frame(ss_output$B_pos[,c(ss_output$f_pos,ss_output$f_neg)])
  plot_df_pos$group = "Stain"
  plot_df_neg = as.data.frame(ss_output$B_neg[,c(ss_output$f_pos,ss_output$f_neg)])
  plot_df_neg$group = "Reference"

  plot_df = rbind(plot_df_pos,plot_df_neg)
  colnames(plot_df) = c("x","y","group")
  ggplot(plot_df, aes(x = x,y = y,col = group)) +
    geom_point()+
    xlab(ss_output$f_pos) + # Label for X-axis
    ylab(ss_output$f_neg) + # Label for Y-axis
    theme_bw()

}
