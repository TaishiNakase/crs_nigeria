#################################################################################
# `post_pred_dist_plot` takes as arguments a vector of states `these_states`, 
#  the output of the MCMC fitting procedure and the raw rubella serology data and 
#  returns a gpplot object with the posterior predictive distribution. 
#################################################################################

post_pred_dist_plot <- function(these_states, MCMC_fit, sero_data, 
                                ng0_MCMC_fit=NULL, ng0_sero_data=NULL, 
                                states) {
  compute_sero_pos_ests <- function(this_state, this_sero_data, this_MCMC_fit) {
    # raw data
    state_sero_data <- this_sero_data %>% 
      filter(state==this_state) %>% 
      mutate(age=age, seropos=positive/tested)
    total_tested <- state_sero_data %>% pull(tested) %>% sum()
    
    # extract posterior predictive distribution
    if (this_state=="Nigeria") state_id <- 1
    else state_id <- which(pull(arrange(filter(this_sero_data, state %in% states), -desc(state)), state)==this_state)[1]
    num_obs <- nrow(state_sero_data)
    pred_pos <- rstan::extract(this_MCMC_fit)[["pos_rep"]][, state_id:(state_id+num_obs-1)]
    
    # summarize the posterior predictive distribution
    num_tested <- this_sero_data %>% filter(state==this_state) %>% arrange(-desc(age)) %>% pull(tested)
    sero_pos_ests <- sapply(1:ncol(pred_pos), function(x) 
    {c(mean(pred_pos[, x]), quantile(pred_pos[, x], c(0.025, 0.975)))/num_tested[x]}) %>%
      t() %>%
      as.data.frame() %>%
      cbind(state_sero_data %>% arrange(-desc(age)) %>% pull(age), .) %>%
      setNames(c("age", "md", "l95", "u95")) %>% 
      mutate(state=this_state)
    return(sero_pos_ests)
  }
  if (all(these_states=="Nigeria")) focal_states <- "Nigeria"
  else focal_states <- setdiff(these_states, "Nigeria")
  sero_pos_ests_list <- purrr::map(focal_states, function(x) 
    compute_sero_pos_ests(this_state=x, this_sero_data=sero_data, this_MCMC_fit=MCMC_fit))
  names(sero_pos_ests_list) <- focal_states
  
  if (("Nigeria" %in% these_states) && !is.null(ng0_MCMC_fit) && !is.null(ng0_sero_data)) {
    sero_pos_ests_list[["Nigeria"]] <- compute_sero_pos_ests(this_state="Nigeria", 
                                                             this_sero_data=ng0_sero_data,
                                                             this_MCMC_fit=ng0_MCMC_fit)
  }
  sero_pos_ests <- do.call(rbind, sero_pos_ests_list) %>% 
    as.data.frame() %>% 
    mutate(state=factor(state, levels=these_states))
  sero_data %<>% 
    rbind(., ng0_sero_data) %>% 
    filter(state %in% these_states) %>% 
    mutate(seropos=positive/tested) %>% 
    mutate(state=factor(state, levels=these_states))

  
  # plot 
  state_label_data <- sero_data %>% 
    group_by(state) %>%  
    summarise(total_tested=sum(tested)) %>% 
    as.data.frame()
  state_labels <- as.character(pull(state_label_data, state))
  state_labels[which(state_labels=="Federal Capital Territory")] <- "FCT"
  state_labels <- paste0(state_labels, " (N=", 
                         formatC(pull(state_label_data, total_tested), big.mark=","),
                         ")")
  names(state_labels) <- pull(state_label_data, state)
  p <- ggplot() + 
    geom_point(data=sero_data, aes(x=age, y=seropos), size=0.5) + 
    geom_line(data=sero_pos_ests, aes(x=age, y=md), color="steelblue3") + 
    geom_ribbon(data=sero_pos_ests,
                aes(ymin=l95, ymax=u95, x=age), fill="steelblue3", alpha = 0.3) + 
    scale_x_continuous(limits=c(0, 45), breaks=seq(0, 45, 5)) + 
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, 0.2)) + 
    labs(x="Age (years)", y="Seroprevalence") + 
    theme_bw() + 
    theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
          panel.background=element_blank(), panel.border=element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(legend.position="none") + 
    theme(axis.text = element_text(size=9), axis.title.y = element_text(size=9), 
          axis.title.x=element_blank(), strip.text=element_text(size=9), 
          strip.background=element_blank()) + 
    facet_rep_wrap(~state, ncol=4, labeller=labeller(state=state_labels))
    
  return(p)
}