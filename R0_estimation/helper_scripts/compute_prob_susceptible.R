#################################################################################
# `compute_prob_sus` computes the probability that an individual at age a years
# is susceptible to infection (i.e. no past history of rubella infection). It 
# takes as arguments the age of the individual in years `a`, the age-specific 
# force of infection estimates `lambdas`, the age groups `agegrps` and the 
# estimate for the rate of loss of maternally-derived maternal immunity. 
#################################################################################

# compute the probability that an individual is susceptible at age a
compute_prob_susceptible <- function(a, lambdas, agegrps, delta) {

  # computes the first expression in the equation
  compute_expression1 <- function(a, lambdas, agegrps, delta) {
    left_endpoint_id <- max(which(agegrps<=a))
    if (left_endpoint_id==1) exp1 <- (1/(lambdas[1]-delta))*exp((lambdas[1]-delta)*a)
    else {
      tmp_sum <- 0 
      for (s in 1:(left_endpoint_id-1)) {
        tmp_sum <- tmp_sum + (1/(lambdas[s]-delta))*
          exp(compute_cum_hazard(agegrps[s+1], lambdas, agegrps)-delta*agegrps[s+1])
      }
      tmp_sum <- tmp_sum + (1/(lambdas[left_endpoint_id]-delta))*
        exp(compute_cum_hazard(agegrps[left_endpoint_id], lambdas, agegrps)+
              lambdas[left_endpoint_id]*(a-agegrps[left_endpoint_id])-delta*a)
      exp1 <- tmp_sum
    }
    return(exp1)
  }
  
  # compute the second expression in the equation
  compute_expression2 <- function(a, lambdas, agegrps, delta) {
    left_endpoint_id <- max(which(agegrps<=a))
    if (left_endpoint_id==1) exp2 <- 0
    else {
      exp2 <- 0 
      for (s in 2:left_endpoint_id) {
        exp2 <- exp2 + (1/(lambdas[s]-delta))*
          exp(compute_cum_hazard(a=agegrps[s], lambdas, agegrps)-delta*agegrps[s])
      }
    }
    return(exp2)
  }
  
  # compute probability susceptible
  prob_sus <- (delta/(delta-lambdas[1]))*exp(-compute_cum_hazard(a=a, lambdas=lambdas, agegrps=agegrps))*
    (1+(delta-lambdas[1])*(compute_expression1(a=a, lambdas=lambdas, agegrps=agegrps, delta=delta)-
                             compute_expression2(a=a, lambdas=lambdas, agegrps=agegrps, delta=delta)))
  return(prob_sus)
}