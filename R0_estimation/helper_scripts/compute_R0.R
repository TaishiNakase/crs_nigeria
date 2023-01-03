#################################################################################
# Functions required for the estimation of the basic reproductive number. 
#################################################################################


#################################################################################
# `compute_I` computes the number of infectious individual in each age group 
#  at equilibrium. It takes as arguments the force of infection vector `lambdas`, 
#  the mortality rate vector `mus`, the rate of loss of maternally-derived 
#  passive immunity, the integer index of the age group `agegrp`, 
#  the vector of age groups `agegrps`, the initial number of individuals `initN`, 
#  the total number of individuals `N`, the number of individuals at equilibrium
#  `equN` and the average duration of infection `D`. 
#################################################################################

compute_I <- function(lambdas, mus, delta, agegrp, agegrps, initN, N, equN, D) {
  # compute life expectancy
  compute_alive_prob <- function(mus, agegrps, a) {
    get_mu <- function(x) {
      return(mus[max(which(agegrps<=x))])
    }
    vec_get_mu <- Vectorize(FUN=get_mu, vectorize.args="x")
    return(exp(-integrate(vec_get_mu, lower=0, upper=a)$value))
  }
  vec_compute_alive_prob <- Vectorize(FUN=compute_alive_prob, vectorize.args="a")
  L <- integrate(vec_compute_alive_prob, mus=mus, agegrps=agegrps, lower=0, upper=Inf)$value
  
  compute_prob_sus_star <- function(x) {
    return(compute_prob_susceptible(a=x, lambdas=lambdas, agegrps=agegrps, delta=delta)*
             compute_alive_prob(mus=mus, agegrps=agegrps, a=x))
  }
  vec_compute_prob_sus_star <- Vectorize(FUN=compute_prob_sus_star, vectorize.args="x")
  left_ep <- agegrps[agegrp]
  if (agegrp==length(agegrps)) {
    right_ep <- Inf
  } else right_ep <- agegrps[agegrp+1]
  equilibriumS <- N*initN[agegrp]/equN[agegrp]/L*
    integrate(vec_compute_prob_sus_star, lower=left_ep, upper=right_ep)$value
  equilibriumI <- lambdas[agegrp]*D*equilibriumS
  return(equilibriumI)
}

#################################################################################
# `compute_betas` computes the estimates of the effective contact rate. 
#  It takes as arguments the force of infection vector `lambdas`, 
#  the mortality rate vector `mus`, the rate of loss of maternally-derived 
#  passive immunity, the life expectancy `L`, the structure of the contact matrix 
#  `contact_mat`,  the vector of age groups `agegrps`, the initial number of 
#  individuals `initN`, the total number of individuals `N`, the number of 
#  individuals at equilibrium equN` and the average duration of infection `D`. 
#################################################################################

compute_betas <- function(lambdas, mus, delta, agegrps, initN, L, D, N=N, equN=equN, contact_mat) {
  equilibriumI <- sapply(1:(length(agegrps)-1), 
                         function(x) compute_I(lambdas=lambdas, mus=mus, 
                                               delta=delta, agegrp=x, 
                                               agegrps=agegrps[-length(agegrps)], 
                                               initN=initN, N=N, equN=equN, D=D))
  
  D_lambda_mat <- matrix(0, nrow=length(agegrps)-1, ncol=length(agegrps)-1)
  for (ii in 1:nrow(contact_mat)) {
    for (jj in 1:ncol(contact_mat)) {
      D_lambda_mat[ii, contact_mat[ii, jj]] <- D_lambda_mat[ii, contact_mat[ii, jj]] + equilibriumI[jj]
    }
  }
  inv_D_lambda_mat <- solve(D_lambda_mat)
  
  betas <- (inv_D_lambda_mat %*% lambdas) %>% as.numeric()
  betas[which(betas<0)] <- 0 # non-regular configurations of beta matrix
  
  return(betas)
}

#################################################################################
# `compute_M_mat` computes the M matrix neccessary for estimation of R0. 
# It takes as arguments the mortality rate vector `mu`, the structure of the 
#  contact matrix `contact_mat` and the vector of age groups `agegrps`,
#################################################################################

compute_M_mat <- function(mu, contact_mat, agegrps) {
  M_mat <- matrix(0, nrow=nrow(contact_mat), ncol=ncol(contact_mat))
  for (ii in 1:nrow(contact_mat)) {
    for (jj in 1:ncol(contact_mat)) {
      if (ii==jj) {
        M_mat[ii, jj] <- compute_avg_life_exp(mu=mu, agegrps=agegrps[-length(agegrps)], 
                                              interval=c(agegrps[ii], agegrps[ii+1]))
      }
    }
  }
  return(M_mat)
}

#################################################################################
# `compute_R0` computes the basic reproductive number. It takes as arguments 
# the vector of age groups `agegrps`, the life expectancy `L`, the structure of 
# the contact matrix `contact_mat`,  the vector of age groups `agegrps`, 
# the initial number of individuals `initN`, the total number of individuals `N`, 
# the number of ndividuals at equilibrium equN`, the average duration of 
# infection `D`, the estimates of effective contact rates `betas` and the 
#  M matrix `M_mat`. 
#################################################################################

compute_R0 <- function(agegrps, initN, N, equN, L, D, contact_mat, betas, M_mat) {
  beta_mat <- matrix(0, nrow=nrow(contact_mat), ncol=ncol(contact_mat))
  for (ii in 1:nrow(contact_mat)) {
    for (jj in 1:ncol(contact_mat)) {
      beta_mat[ii, jj] <- betas[contact_mat[ii, jj]]
    }
  }
  B_mat <- t(apply((N*D/L)*beta_mat, 1, function(x) x*initN/equN))
  NGM <- M_mat %*% B_mat
  R0 <- max(eigen(NGM, only.values=TRUE)$values)
  return(R0)
}