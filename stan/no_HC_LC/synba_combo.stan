functions{

  real synba(real d1, real d2, real E0, real E1, real E2, real E3, real C1, real C2, real h1, real h2, real a){
    real d1h1;
    real d2h2;
    real C1h1;
    real C2h2;
    real u;

    C1h1 = pow(C1,h1);
    C2h2 = pow(C2,h2);
    d1h1 = pow(d1,h1);
    d2h2 = pow(d2,h2);
    u = (C1h1*C2h2*E0 + d1h1*C2h2*E1 + C1h1*d2h2*E2 + a*d1h1*d2h2*E3) / (C1h1*C2h2 + d1h1*C2h2 + C1h1*d2h2 + a*d1h1*d2h2);
//    Y[i] ~ normal(((exp(logC_1))^h_1 * (exp(logC_2))^h_2 * e_0 + X1[i]^h_1 * (exp(logC_2))^h_2 * e_1 * e_0 + X2[i]^h_2 * (exp(logC_1))^h_1 * e_2 * e_0 + alpha * X1[i]^h_1 * X2[i]^h_2 * e_3 * e_0) / ((exp(logC_1))^h_1 * (exp(logC_2))^h_2 + X1[i]^h_1 * (exp(logC_2))^h_2 + X2[i]^h_2 * (exp(logC_1))^h_1 + alpha * X1[i]^h_1 * X2[i]^h_2), sigma);

    return( u );
  }

}


data{

  // Observed variables
  int<lower=0> Nc; // The number of data points for combination drug
  array[Nc] real yc; // target cell viability with combination drug
  array[Nc] real<lower=0> dc1; // Concentrations for drug1 in combination
  array[Nc] real<lower=0> dc2; // Concentrations for drug2 in combination

  // Fixed parameters
  real<lower=0> hc_mean; // Target cell viability without treatment
  real<lower=0> einf_beta_a;
  real<lower=0> einf_beta_b;
  real lb_log10_C;
  real ub_log10_C;
  real sigma_mu;

  // Variables for prediction
  int<lower=0> Nc_new; // The number of data points
  array[Nc_new] real<lower=0> dc1_new; // Concentrations
  array[Nc_new] real<lower=0> dc2_new; // Concentrations

}


parameters{

  // For basic dose response curve
  real<lower=0> E0;
  real<lower=0, upper=1> e_1;
  real<lower=0, upper=1> e_2;
  real<lower=0, upper=1> e_3;
  real<lower=lb_log10_C, upper=ub_log10_C> log_C1;
  real<lower=lb_log10_C, upper=ub_log10_C> log_C2;
  real<lower=0> h1;
  real<lower=0> h2;
  real<lower=0> s_y;

  // For potency
  real<lower=0> alpha;
}


transformed parameters{

  real<lower=0> E1;
  real<lower=0> E2;
  real<lower=0> E3;
  real<lower=0> C1;
  real<lower=0> C2;

  // Efficacy
  real beta1;
  real beta2;

  E1 = E0 * e_1;
  E2 = E0 * e_2;
  E3 = E0 * e_3;

  C1 = pow(10, log_C1);
  C2 = pow(10, log_C2);

  beta1 = (E1 - E3)/(E0 - E1);
  beta2 = (E2 - E3)/(E0 - E2);
}


model{

  // Expected target cell viability
  array[Nc] real uc;

  E0 ~ normal(hc_mean, 0.03*hc_mean);

  e_1 ~ beta(einf_beta_a, einf_beta_b);
  e_2 ~ beta(einf_beta_a, einf_beta_b);
  e_3 ~ beta(einf_beta_a, einf_beta_b);

  alpha ~ lognormal(0, 1);
  h1 ~ lognormal(0, 1);
  h2 ~ lognormal(0, 1);
  s_y ~ lognormal(sigma_mu, 1);

  for (n in 1:Nc) {
    uc[n] = synba(dc1[n], dc2[n], E0, E1, E2, E3, C1, C2, h1, h2, alpha);
    target += normal_lpdf(yc[n] | uc[n], s_y);
  }

}


generated quantities{

  // Expected target cell viability
  array[Nc_new] real<lower=0> uc_pred;
  real<lower=0> u_max;

  // Predicted target cell viability
  array[Nc_new] real yc_pred;
  array[Nc] real yc_hat;

  // Pointwise log likelihood
  vector[Nc] log_lik;

  u_max = E0;

  for (n in 1:Nc_new) {
    uc_pred[n] = synba(dc1_new[n], dc2_new[n], E0, E1, E2, E3, C1, C2, h1, h2, alpha);
    yc_pred[n] = normal_rng(uc_pred[n], s_y);
  }

  for (n in 1:Nc) {
    yc_hat[n] = normal_rng(synba(dc1[n], dc2[n], E0, E1, E2, E3, C1, C2, h1, h2, alpha), s_y);
    log_lik[n] = normal_lpdf(yc[n] | synba(dc1[n], dc2[n], E0, E1, E2, E3, C1, C2, h1, h2, alpha), s_y);
  }

}

