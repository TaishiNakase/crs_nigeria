#################################################################################
# `plot_foi_distribution` takes as arguments a state name `this_state`, 
#  a matrix of lambda values and the vector of the age groups, 
#  and returns a ggplot object for the force of infection across age 
#  with 90% credible intervals.
#################################################################################

foi_distribution_plot <- function(this_state, lambdas, age_groups, ymax=1) {
  f <- function(xmin, xmax, lambda) {
    return(data.frame(age=seq(xmin, xmax, 0.1), mean=mean(lambda), 
                      l95=as.numeric(quantile(lambda, 0.025)), u95=as.numeric(quantile(lambda, 0.975))))
  }
  df <- do.call(rbind, purrr::map(1:(length(age_groups)-1), function(x) 
    f(xmin=age_groups[x], xmax=age_groups[x+1], 
      lambda=lambdas[, x]))) %>%
    filter(age<=45)
  
  p <- ggplot(df) + 
    geom_line(aes(x=age, y=mean), color="steelblue3") + 
    geom_ribbon(aes(ymin=l95, ymax=u95, x=age), fill="steelblue3", alpha = 0.3) + 
    scale_x_continuous(limits=c(0, 45), breaks=seq(0, 45, 5)) + 
    scale_y_continuous(limits=c(0, ymax), breaks=seq(0, ymax, 0.1)) + 
    labs(x="Age (year)", y="Force of infection (1/yr)", title=this_state) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
    theme(legend.position="none")
  return(p)
}