// user-defined functions
functions {
  // computes the cumulative hazard for a given age
  real compute_cum_hazard(real a, real[] lambdas, real[] agegrps) {
    int left_endpoint_id;
    real left_endpoint_agegrp;
    real cum_hazard = 0;
    
    // identify age group 
    for (i in 1:num_elements(agegrps)) {
      if (a>=agegrps[i]) left_endpoint_id = i;
      else break;
    }
    left_endpoint_agegrp = agegrps[left_endpoint_id];
    
    // compute cumulative hazard
    if (left_endpoint_id==1) cum_hazard = lambdas[1]*a;
    else {
      for (i in 1:(left_endpoint_id-1)) {
        cum_hazard += (agegrps[i+1]-agegrps[i])*lambdas[i];
      }
      cum_hazard += lambdas[left_endpoint_id]*(a-left_endpoint_agegrp);
    }
    return(cum_hazard);
  }
  
  // computes the first expression of the probability of susceptible as a function of age 
  real compute_expression1(real a, real[] lambdas, real[] agegrps, real delta) {
    int left_endpoint_id;
    real exp1 = 0;
    
    // identify age group 
    for (i in 1:num_elements(agegrps)) {
      if (a>=agegrps[i]) left_endpoint_id = i;
      else break;
    }
    
    // compute first expression
    if (left_endpoint_id==1) exp1 = (1/(lambdas[1]-delta))*exp((lambdas[1]-delta)*a);
    else {
      for (s in 1:(left_endpoint_id-1)) {
        exp1 += (1/(lambdas[s]-delta))*exp(compute_cum_hazard(agegrps[s+1], lambdas, agegrps)-delta*agegrps[s+1]);
      }
      exp1 += (1/(lambdas[left_endpoint_id]-delta))*exp(compute_cum_hazard(agegrps[left_endpoint_id], lambdas, agegrps)+lambdas[left_endpoint_id]*(a-agegrps[left_endpoint_id])-delta*a);
    }
    return(exp1);
  }
  
  // computes the second expression of the probability of susceptible as a function of age 
  real compute_expression2(real a, real[] lambdas, real[] agegrps, real delta) {
    int left_endpoint_id;
    real express2 = 0;
    
    for (i in 1:num_elements(agegrps)) {
      if (a>=agegrps[i]) left_endpoint_id = i;
      else break;
    }
    
    if (left_endpoint_id==1) express2 = 0;
    else {
      for (s in 2:left_endpoint_id) {
        express2 += (1/(lambdas[s]-delta))*exp(compute_cum_hazard(agegrps[s], lambdas, agegrps)-delta*agegrps[s]);
      }
    }
    return(express2);
  }
}

// data variables 
data {
  int<lower=0> N;          // total num measurements
  int<lower=0> K;          // number of states
  int<lower=0> gs[K];      // number of observations for each state vector
  int<lower=0> pos[N];     // num seropositive vector
  int<lower=0> tested[N];  // num tested vector
  real<lower=0.0> age[N];  // age vector
  int<lower=0> num_agegrp; // num age groups
  int agegrp_id[N];        // age group id vector
  real agegrp[num_agegrp]; // age groups vector
  int group_id[K];         // group id vector
  int<lower=0> num_groups; // number of unique groups
}

// transformed data 
transformed data {
  int<lower=0> st_pos[K];  // start position for state data in vectors
  for (i in 1:K) {
    if (i==1) st_pos[i] = 1;
    else st_pos[i] = st_pos[i-1] + gs[i-1];
  }
}

// the parameters
parameters {
  real<lower=0.0> lambda[num_agegrp*K-(K-num_groups)];
  real<lower=0.0> delta_inv;
}

transformed parameters {
  real<lower=0.0, upper=1.0> p[N];
  real<lower=0.0> delta = 1./delta_inv;
  for (i in 1:K) {
    real state_lambdas[num_agegrp]; 
    for (k in 1:(num_agegrp-1)) {
      state_lambdas[k] = lambda[K*(k-1)+i];
    }
    state_lambdas[num_agegrp] = lambda[K*(num_agegrp-1)+group_id[i]];
    
    for (j in 1:gs[i]) {
      int s = st_pos[i]+(j - 1);
      p[s] = 1-(delta/(delta-state_lambdas[1]))*exp(-compute_cum_hazard(age[s], state_lambdas, agegrp))*
             (1+(delta-state_lambdas[1])*(compute_expression1(age[s], state_lambdas, agegrp, delta)-
             compute_expression2(age[s], state_lambdas, agegrp, delta)));
    }
  }
}

// model
model {
  lambda ~ lognormal(-2.6, 0.7);
  delta_inv ~ lognormal(-1.04, 0.74);
  for (i in 1:N) {
    pos[i] ~ binomial(tested[i], p[i]);
  }
}

// generated quantities
generated quantities {
  real pos_rep[N];
  pos_rep = binomial_rng(tested, p);
}
