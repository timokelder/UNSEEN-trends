#' Internal function to plot histograms for the boottest
#' @param bootstrapped_ens
#' @param bootstrapped_ensratio
#' @param obs
#' @param fun
#' @param main
#' @param units
#' @param fontsize
#' @noRd # no documentation
plot_hist_combined <- function(bootstrapped_ens, bootstrapped_ensratio, obs, fun, main, units, fontsize) {
  bootstrapped_fun <- apply(bootstrapped_ens, MARGIN = 2, FUN = fun)
  bootstrapped_fun_ratio <- apply(bootstrapped_ensratio, MARGIN = 2, FUN = fun)

  ggplot2::ggplot() +
    ggplot2::geom_histogram(ggplot2::aes(x = bootstrapped_fun), 
                            color = "black", fill = "black", 
                            alpha = 0.5,
                            bins = 30) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(bootstrapped_fun, 
                                                           probs = c(0.025, 0.975))),
                        color = "black", linetype = "dashed", size = 1) +
    ggplot2::geom_histogram(ggplot2::aes(x = bootstrapped_fun_ratio),
                            color = "black", fill = "orange",
                            alpha = 0.5, bins = 30) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = quantile(bootstrapped_fun_ratio,
                                         probs = c(0.025, 0.975))),
                        color = "orange", linetype = "dashed", size = 1) +
    ggplot2::geom_vline(ggplot2::aes(xintercept = fun(obs)),
                        color = "blue", size = 2
    ) +
    ggplot2::labs(title = main, y = "Number of bootstrapped series", x = paste0(" (", units, ")")) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      text = ggplot2::element_text(size = fontsize),
      axis.text = ggplot2::element_text(size = fontsize)
    )
}


#' Test the model stability: 1. density plot
#'
#' We plot the density distribution of the different leadtimes using ggplot.
#'
#'
#' @param obs     The observations. This function expects a vector, i.e. dataframe$variable. 
#' @param ensemble The UNSEEN ensemble. This function expects a vector, i.e. dataframe$variable
#' @param fontsize 
#'
#' @return plots showing the bootstrapped tests of the mean, sd, skewness and kurtosis
#' @source Evaluation explaned in more detail in Kelder et al. 2020
Boottest_biascor <- function(obs, ensemble, fontsize = 11) {
  bootstrapped_ens <- sample(ensemble, size = length(obs) * 10000, replace = T) # The original raw data
  bootstrapped_ens <- array(bootstrapped_ens, dim = c(length(obs), 10000)) # Creates an array with 10.000 series of 35 values
  ensemble_ratio <- ensemble * mean(obs) / mean(ensemble) ## The simple ratio as mean bias corrected series
  bootstrapped_ensratio <- base::sample(ensemble_ratio, size = length(obs) * 10000, replace = T) # The original raw data
  bootstrapped_ensratio <- array(bootstrapped_ensratio, dim = c(length(obs), 10000)) # Creates an array with 10.000 series of 35 values

  p1_comb <- plot_hist_combined(bootstrapped_ens, bootstrapped_ensratio, obs, fun = mean, main = "", units = "mm/day", fontsize = fontsize) # Mean
  p2_comb <- plot_hist_combined(bootstrapped_ens, bootstrapped_ensratio, obs, fun = sd, main = "", units = "mm/day", fontsize = fontsize) # Standard Deviation
  p3_comb <- plot_hist_combined(bootstrapped_ens, bootstrapped_ensratio, obs, fun = moments::skewness, main = "", units = "-", fontsize = fontsize) # Skewness
  p4_comb <- plot_hist_combined(bootstrapped_ens, bootstrapped_ensratio, obs, fun = moments::kurtosis, main = "", units = "-", fontsize = fontsize) # Kurtosis

  ggpubr::ggarrange(p1_comb, p2_comb, p3_comb, p4_comb,
    labels = c("a", "b", "c", "d"),
    font.label = list(size = fontsize, color = "black", face = "bold", family = NULL),
    ncol = 2, nrow = 2
  )
}


#   p1 <-
#     ggplot2::ggplot(ensemble,
#                     ggplot2::aes(x = precipitation,
#                                  colour = leadtime)) +
#     # ggplot2::ggtitle("UK") +
#     ggplot2::labs(x = lab, y = "Density") +
#     ggplot2::geom_line(stat = "density") +
#     ggplot2::theme_classic() +
#     ggplot2::theme(legend.position = "none") +
#     ggplot2::scale_colour_manual(values = cbPalette) #+
#   # theme(
#   #   text = element_text(size = 11),
#   #   axis.text = element_text(size = 11),
#   #   plot.title = element_text(hjust = 0.5)
#   # )
#
#   return(p1)
# }
