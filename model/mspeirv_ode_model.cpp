#include <numeric>
#include <stdio.h>
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector conduct_sia(NumericVector& y, int st_agegrp_idx, int end_agegrp_idx, double cov) {
  int M_id = y.findName("M1");
  int S_id = y.findName("S1");
  int E_id = y.findName("E1");
  int R_id = y.findName("R1");
  int V_id = y.findName("V1");
  int SP1_id = y.findName("1SP1");
  int SP2_id = y.findName("2SP1");
  int SP3_id = y.findName("3SP1");
  int EP_id = y.findName("EP1");
  
  for (int i = st_agegrp_idx; i <= (end_agegrp_idx); ++i) {
    y[V_id+i-1] = y[V_id+i-1]+cov*(y[M_id+i-1]+y[S_id+(i-1)]+y[SP1_id+(i-1)]+y[SP2_id+(i-1)]+y[SP3_id+(i-1)]+y[E_id+(i-1)]+y[EP_id+(i-1)]+y[E_id+i-1]+y[R_id+i-1]);
    y[M_id+i-1] = (1-cov)*y[M_id+i-1];
    y[S_id+i-1] = (1-cov)*y[S_id+i-1];
    y[SP1_id+i-1] = (1-cov)*y[SP1_id+i-1];
    y[SP2_id+i-1] = (1-cov)*y[SP2_id+i-1];
    y[SP3_id+i-1] = (1-cov)*y[SP3_id+i-1];
    y[E_id+i-1] = (1-cov)*y[E_id+i-1];
    y[EP_id+i-1] = (1-cov)*y[EP_id+i-1];
    y[R_id+i-1] = (1-cov)*y[R_id+i-1];
  }
  return y;
}

// [[Rcpp::export]]
List get_ode_rates(double t, NumericVector& y, NumericVector& parms) {
  // constants
  static double PI = 3.14159265358979323846;
  
  // enum for variable positions
  enum Class {MID, SID, EID, IID, RID, VID, SP1ID, SP2ID, SP3ID, EPID};
  
  // parameters
  static int mu_id = parms.findName("mu1");
  static int b_id = parms.findName("b1");
  static int omega_id = parms.findName("omega1");
  static int sigma_id = parms.findName("sigma1");
  static int theta_id = parms.findName("theta1");
  static int gamma_id = parms.findName("gamma1");
  static int beta_id = parms.findName("beta1");
  static int alpha_id = parms.findName("alpha1");
  static int f_id = parms.findName("f1");
  static int delta1_id = parms.findName("1delta1");
  static int delta2_id = parms.findName("2delta1");
  static int delta3_id = parms.findName("3delta1");
  static int pv_id;
  static int eps_id;
  static int num_stages = parms["num_stages"];
  static int num_classes = parms["num_classes"];
  double n_id = parms.findName("N1");
  int sim_type = parms["sim_type"];
  double db = parms["db"];
  double mu, b, omega, sigma, theta, gamma, f, delta1, delta2, delta3, pv=0.0, eps, initN;
  static NumericVector beta;
  static NumericVector alpha;
  double seasonal_force = cos(2*PI*t/365);
  
  // variables
  static int M_id = y.findName("M1");
  static int S_id = y.findName("S1");
  static int E_id = y.findName("E1");
  static int I_id = y.findName("I1");
  static int R_id = y.findName("R1");
  static int V_id = y.findName("V1");
  static int SP1_id = y.findName("1SP1");
  static int SP2_id = y.findName("2SP1");
  static int SP3_id = y.findName("3SP1");
  static int EP_id = y.findName("EP1");
  static double M, S, E, I, R, V, SP1, SP2, SP3, EP, aN;
  double N = sum(y[seq(0, num_classes*num_stages-1)]);
  NumericVector all_I = y[seq(0, num_stages-1)+IID*num_stages];
  
  // intermediates
  static NumericVector rates (num_stages*num_classes);
  static double FOI;
  int yr = (int)(t/365);
  static double dM, dS, dE, dI, dR, dV, dSP1, dSP2, dSP3, dEP;
  
  // routine vaccination
  if (sim_type!=0) pv_id = parms.findName("pv1");
  
  // external infected immigration
  eps_id = parms.findName("eps");
  eps = parms[eps_id];
  
  for (int i = 0; i < num_stages; ++i) {
    // variables
    M = y[M_id+i];
    S = y[S_id+i];
    E = y[E_id+i];
    I = y[I_id+i];
    R = y[R_id+i];
    V = y[V_id+i];
    SP1 = y[SP1_id+i];
    SP2 = y[SP2_id+i];
    SP3 = y[SP3_id+i];
    EP = y[EP_id+i];
    aN = M+S+E+I+R+V+SP1+SP2+SP3+EP;
    
    // parameters
    initN = parms[n_id+i];
    mu = parms[mu_id+i];
    b = parms[b_id+i]*pow(1+db, t);
    omega = parms[omega_id+i];
    sigma = parms[sigma_id+i];
    gamma = parms[gamma_id+i];
    theta = parms[theta_id+i];
    f = parms[f_id+i]*pow(1+db, t)*initN/aN;
    delta1 = parms[delta1_id+i];
    delta2 = parms[delta2_id+i];
    delta3 = parms[delta3_id+i];
    beta = parms[seq(0, num_stages-1)*num_stages+i+beta_id];
    alpha = parms[seq(0, num_stages-1)*num_stages+i+alpha_id];
    
    // maternally immune
    dM = b*N - (omega+mu+theta)*M;
    
    // susceptible
    FOI = sum((1+alpha*seasonal_force)*beta*all_I)*initN/aN;
    dS = omega*M + delta3*SP3 - (FOI+mu+theta+f)*S;
    
    // susceptible pregnancies
    dSP1 = f*S - (mu+theta+delta1+FOI)*SP1;
    dSP2 = delta1*SP1 - (mu+theta+delta2+FOI)*SP2;
    dSP3 = delta2*SP2 - (mu+theta+delta3+FOI)*SP3;
    
    // exposed
    dE = (FOI)*(S+SP2+SP3) - (sigma+mu+theta)*E;
    dEP = (FOI)*SP1 - (sigma+mu+theta)*EP;
    
    // infectious
    dI = sigma*(E+EP) + eps*aN - (gamma+mu+theta)*I;
    
    // recovered
    dR = gamma*I - (mu+theta)*R;
    
    // vaccinated
    dV = -(mu+theta)*V;
    
    if (i != 0) {
      theta = parms[theta_id+i-1];
      
      // routine vaccination: starts at 9 mnths (first cohort in 9-10 agegrp)
      if (sim_type!=0 && i==9) {
        pv = parms[pv_id+yr];
      }
      
      dM += (1-pv)*theta*y[M_id+(i-1)];
      dS += (1-pv)*theta*y[S_id+(i-1)];
      dSP1 += (1-pv)*theta*y[SP1_id+(i-1)];
      dSP2 += (1-pv)*theta*y[SP2_id+(i-1)];
      dSP3 += (1-pv)*theta*y[SP3_id+(i-1)];
      dE += (1-pv)*theta*y[E_id+(i-1)];
      dEP += (1-pv)*theta*y[EP_id+(i-1)];
      dI += theta*y[I_id+(i-1)];
      dR += (1-pv)*theta*y[R_id+(i-1)];
      dV += pv*theta*(y[M_id+(i-1)]+y[S_id+(i-1)]+y[SP1_id+(i-1)]+y[SP2_id+(i-1)]+y[SP3_id+(i-1)]+y[E_id+(i-1)]+y[EP_id+(i-1)]+y[R_id+(i-1)])+theta*y[V_id+(i-1)];
      
      pv=0; // routine vaccination
    }
    rates[i+MID*num_stages] = dM;
    rates[i+SID*num_stages] = dS;
    rates[i+EID*num_stages] = dE;
    rates[i+IID*num_stages] = dI;
    rates[i+RID*num_stages] = dR;
    rates[i+VID*num_stages] = dV;
    rates[i+SP1ID*num_stages] = dSP1;
    rates[i+SP2ID*num_stages] = dSP2;
    rates[i+SP3ID*num_stages] = dSP3;
    rates[i+EPID*num_stages] = dEP;
  } 
  return List::create(rates); 
}