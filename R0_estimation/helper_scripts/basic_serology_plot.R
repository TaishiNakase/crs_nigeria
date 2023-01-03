#################################################################################
# `basic serology plot` takes as arguments a state name `this_state` and a 
# dataframe of serological data (age, positive and tested) and returns a plot
# of its serological profile. 
#################################################################################
basic_serology_plot <- function(this_state, sero_data) {
  state_sero_data <- sero_data %>% 
      filter(state==this_state) %>% 
      mutate(seropos=positive/tested)
  total_tested <- state_sero_data %>% pull(tested) %>% sum()

  p <- ggplot(data=state_sero_data) + 
    geom_point(aes(x=age, y=seropos)) + 
    scale_x_continuous(limits=c(0, 45), breaks=seq(0, 45, 5)) + 
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) + 
    labs(x="Age (years)", y="Seroprevalence", 
         title=paste0(this_state, " (N=", total_tested, ")")) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          title = element_text(size=10))
  return(p)
}
  