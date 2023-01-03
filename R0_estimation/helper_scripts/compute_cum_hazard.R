#################################################################################
# `compute_cum_hazard` computes the cummulative hazard of infection for an 
#  `a`-year-old individual. It takes as arguments the age of the individual 
#  in years `a`, the age-specific force of infection estimates `lambdas` and
#  the age groups `agegrps`.
#################################################################################

# computes the cumulative hazard at a given age
compute_cum_hazard <- function(a, lambdas, agegrps) {
  left_endpoint_id <- max(which(agegrps<=a))
  left_endpoint_agegrp <- agegrps[left_endpoint_id]
  if (left_endpoint_id==1) cum_hazard <- lambdas[1]*a
  else {
    cum_hazard <- sum(diff(agegrps[1:left_endpoint_id])*lambdas[1:(left_endpoint_id-1)])
    cum_hazard <- cum_hazard + lambdas[left_endpoint_id]*(a-left_endpoint_agegrp)
  }
  return(cum_hazard)
}