#################################################################################
# `compute_avg_life_exp` computes the average life expectancy. It takes as 
# arguments the vector of age-specific mortality rates, the vector of age groups 
# and an optional interval over which to compute the avg. life exp. 
#################################################################################

compute_avg_life_exp <- function(mu, agegrps, interval=c(0, Inf)) {
  # computes the ccdf of the survival time
  ccdf <- function(x) {
    left_endpoints <- agegrps[(which(agegrps<=x))]
    intervals <- diff(c(left_endpoints, x))
    this_mu <- mu[seq_along(intervals)]
    hazard <- sum(this_mu*intervals)
    return(exp(-hazard))
  }
  vec_ccdf <- Vectorize(FUN=ccdf, vectorize.args="x")
  
  # compute average life expectancy
  avg_life_exp <- integrate(vec_ccdf, lower=interval[1], upper=interval[2])$value
  return(avg_life_exp)
}