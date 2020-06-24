#' Test the model stability: 1. density plot
#' 
#' We plot the density distribution of the different leadtimes using ggplot.   
#'    
#' 
#' @param ensemble The UNSEEN ensemble. This function expects an dataframe with variables leadtime, precipitation.
#'
#' @return a plot showing the empirical probability density distribution for each leadtime 
#' @source Evaluation explaned in more detail in Kelder et al. 2020
#' @source Colorblind friendly palette  http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
Model_stability_density <- function(ensemble, lab = "") {
  # I select five colors and put black at the end.
  cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#D55E00", "#000000") # , "#0072B2", "#CC79A7")
  ### The Leadtime column has to be a factor for a grouped plot
  ensemble$leadtime <- as.factor(ensemble$leadtime)
  
  p1 <-
    ggplot2::ggplot(ensemble, 
                    ggplot2::aes(x = precipitation, 
                                 colour = leadtime)) +
    # ggplot2::ggtitle("UK") +
    ggplot2::labs(x = lab, y = "Density") +
    ggplot2::geom_line(stat = "density") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "none") +
    ggplot2::scale_colour_manual(values = cbPalette) #+
  # theme(
  #   text = element_text(size = 11),
  #   axis.text = element_text(size = 11),
  #   plot.title = element_text(hjust = 0.5)
  # )
  
  return(p1)
}

#' Test the model stability: 2. empirical extreme value plot
#' 
#' We plot the density distribution of the different leadtimes using ggplot.   
#' We want to show the confidence interval of the distribution of all lead times pooled together,
#' and test whether the individual lead times fall within these confidence intervals. 
#' Therefore we bootstrap the pooled leadtimes into series with equal length of the individual leadtimes (875), with n=10000.  
#'    
#' @param ensemble The UNSEEN ensemble. This function expects an dataframe with variables leadtime, precipitation.
#'
#' @return a plot with the empirical return values of the pooled ensemble including confidence intervals. 
#' Individual lead times are plotted on top. 
#' @source Evaluation explaned in more detail in Kelder et al. 2020
#' @source Colorblind friendly palette  http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
Model_stability_boot <- function(ensemble, lab = "") {  
  ## Still to do. See the UNSEEN-open project.
  
  # Leadtime_length <- sum(ensemble$leadtime == 2)
  # bootstrapped_series <- sample(ensemble$precipitation, size = Leadtime_length * 10000, replace = T) # bootstraps the series of length equal to each lead time (875) with n= 10.000
  # bootstrapped_array <- array(bootstrapped_series, dim = c(Leadtime_length, 10000)) # Creates an array with 10.000 series of 875 values
  # 
  # CI <- function(x) {
  #   quantile(x, probs = c(0.025, 0.975)) ## lower and upper interval
  # 
  
}

