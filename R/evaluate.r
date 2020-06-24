#' Illustrate the timeseries of the UNsEEN ensemble and the observed using ggplot
#'
#' @param ensemble The UNSEEN ensemble. This must be a dataframe with at least the variables year and precipitation
#' @param obs The observations. This must be a dataframe with the variables year and precipitation
#' @param ylab The y label (string)
#' @param title The title (string)
#' @return a ggplot with x = year and y =
#'
plot_unseen <- function(ensemble, obs, ylab = "", title = "") {
  year <- precipitation <- NULL
  ggplot2::ggplot() +
    ggplot2::geom_boxplot(
      data = ensemble,
      mapping = ggplot2::aes(
        x = year, y = precipitation,
        group = year,
        fill = "UNSEEN"
      ),
      alpha = 0.3
    ) + ## Seas5 color is defined manually
    ggplot2::geom_point(
      data = obs,
      ggplot2::aes(
        x = year,
        y = precipitation,
        col = "OBS"
      ),
      shape = 4,
      size = 2,
      stroke = 1.5
    ) +
    ggplot2::theme_classic() +
    ggplot2::scale_fill_manual(name = NULL, values = c("UNSEEN" = "black")) + ## Here SEAS5 color is defined
    ggplot2::scale_colour_manual(name = NULL, values = c("OBS" = "blue")) + ## And ERA5 color
    ggplot2::ylab(ylab) +
    ggplot2::ggtitle(title) +
    ggplot2::theme(
      legend.position = c(.95, .02),
      legend.justification = c("right", "bottom"),
      legend.spacing.y = ggplot2::unit(-0.2, "cm"),
      legend.title = ggplot2::element_blank()
    ) + # ,
    # text=element_text(size=7),
    #   axis.text = element_text(size=7))+
    ggplot2::guides(
      color = ggplot2::guide_legend(order = 1),
      fill = ggplot2::guide_legend(order = 2)
    )
}


#' Calculate the correlation between each unique pair of ensemble members.
#'
#' @param ensemble The UNSEEN ensemble. This function expects an dataframe with variables ensemble, leadtime, precipitation.
#' It expects that the first member is 0 (as used in CDS) and that the first lead time is 2 
#' (first month shows dependence and is therefore removed).
#'
#' @return a plot showing the boxplot of the pairwise correlations for each leadtime. See Kelder et al. 2020
#' @section Warning:
#' The confidence intervals are based on a bootstrap test with n_ensembles = 25. This will not give the right CI bounds for different
#' ensemble sizes. 
Independence_test <- function(ensemble, n_ensembles = 25, n_lds = 5) {
  `2.5%` <- `97.5%` <- Boxstat <- Freq <- cor <- leadtime <- x <- precipitation <- NULL
  

  ## Perform the pairwise test
  correlations_lds <- array( # create and array to fill the correlations in
    dim = c(n_ensembles, n_ensembles, n_lds), # Dimensions
    dimnames = list(
      memberx = as.character(0:(n_ensembles - 1)), # The names of the members start at 0
      membery = as.character(0:(n_ensembles - 1)), # Lead times start at 2 (first month shows dependence and is therefore removed)
      leadtime = as.character(2:(n_lds + 1))
    )
  )
  # lds, mbmrs
  for (ld in 2:(n_lds + 1)) {
    for (mbr1 in 0:(n_ensembles - 1)) {
      for (mbr2 in 0:(n_ensembles - 1)) {
        if (mbr1 > mbr2) { ## Only calculate this for the top half of the correlation matrix, as not to duplicate any values 
                           ## -> avoid correlating mbr 1 with mbr2 and then mbr2 with mbr 1.
          predictant <- ensemble$precipitation[ensemble$ensemble == mbr1 & ensemble$leadtime == ld]
          predictor <- ensemble$precipitation[ensemble$ensemble == mbr2 & ensemble$leadtime == ld]

          correlations_lds[as.character(mbr1), as.character(mbr2), as.character(ld)] <- cor(
            x = predictant,
            y = predictor,
            method = "spearman"
          )
        }
      }
    }
  }

  ##Convert to dataframe for ggplot
  correlations_lds_df = as.data.frame(as.table(correlations_lds)) 
  ## And plot
  CI_bounds_boxplots_df = as.data.frame(t(CI_bounds_boxplots))
  CI_bounds_boxplots_df$Boxstat = c('lower whisker', 'lower IQR', 'median','upper IQR', 'upper whisker')
  CI_bounds_boxplots_df_ld2 <- CI_bounds_boxplots_df
  CI_bounds_boxplots_df_ld2$x <- 0.5
  CI_bounds_boxplots_df_ld_end <- CI_bounds_boxplots_df
  CI_bounds_boxplots_df_ld_end$x <- n_lds + .5
  CI_bounds_boxplots_df_all = rbind(CI_bounds_boxplots_df_ld2,CI_bounds_boxplots_df_ld_end)
  
  ggplot2::ggplot(data = correlations_lds_df) +
    ggplot2::geom_violin(
      mapping = ggplot2::aes(
        x = leadtime, 
        y = Freq,
        group = leadtime
      )
    ) +
    ggplot2::geom_boxplot(mapping = ggplot2::aes(x = leadtime, y = Freq,group = leadtime), width=0.2) + # same as violin
    ggplot2::geom_ribbon(data=CI_bounds_boxplots_df_all, 
                         ggplot2::aes(x = x, 
                                      ymin = `2.5%`,ymax=`97.5%`, 
                                      group = Boxstat),
                         fill="grey",
                         alpha="0.4") +
    ggplot2::theme_classic() +
    ggplot2::ylab(bquote("Spearman" ~ rho))
  
  
# The old plot with base R
#   # par(mfrow=c(1,2),mar=c(4.5,5.1,2.1,2.1),cex.axis=1.5, cex.lab=1.5,cex.main=1.5)
#   boxplot(list(correlations_lds[, , "2"], correlations_lds[, , "3"], correlations_lds[, , "4"], correlations_lds[, , "5"], correlations_lds[, , "6"]),
#     range = 0, # box whiskers to be the data extremes
#     xaxt = "n",
#     xlab = "Lead time",
#     ylab = bquote("Spearman" ~ rho)
#   )
#   for (i in 1:length(CI_bounds_boxplots[1, ])) {  ## CI_bounds_boxplots is internal data, created by function b
#     polygon(c(0, 6, 6, 0), c(rep(CI_bounds_boxplots[1, i], 2), rep(CI_bounds_boxplots[2, i], 2)), col = gray(0.8, alpha = 0.3))
#   }
#   Axis(side = 1, at = 1:6, labels = c(as.character(2:6), "all"))

    # return(correlations_lds)
}


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

