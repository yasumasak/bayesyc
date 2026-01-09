
functions{
  real logistic4(real d, real E0, real E1, real C, real h){
    return( E1 + (E0 - E1) / (1 + pow(d/C, h)) );
  }
}

data{
  // Observed variables
  int<lower=0> N; // The number of data points
  array[N] int<lower=0> y; // target cell count
  array[N] real<lower=0> x; // Concentrations

  int<lower=0> N_hc; // The number of high control data points
  array[N_hc] int<lower=0> hc; // High control target cell count
  real<lower=0> hc_med; // Median of high control target cell count

  int<lower=0> N_lc; // The number of low control data points
  array[N_lc] int<lower=0> lc; // Low control target cell count

  // hyperparameters
  real lb_log10_C; // Lower bound for log10(C)
  real ub_log10_C; // Upper bound for log10(C)

  // Variables for prediction
  int<lower=0> N_new; // The number of data points
  array[N_new] real<lower=0> x_new; // Concentrations
}

parameters{
  // For dose response curve
  real<lower=0> E0;
  real<lower=0, upper=1> e1;
  real log10_C;
  real<lower=0> h;
  real<lower=0> noise;
  real<lower=0> s_lc;
  real<lower=0> s_y;

  // For outlier
  array[N] real<lower=0, upper=1> pi; // normal: 0, outlier: 1
  array[N_hc] real<lower=0, upper=1> pi_hc; // normal: 0, outlier: 1
  real<lower=0> k;
}

transformed parameters{
  real<lower=0> Einf;
  real<lower=0> C;
  array[N] real<lower=0> u;

  // For E0 > Einf
  Einf = E0 * e1;
  C = pow(10, log10_C);

  for (n in 1:N) {
    u[n] = logistic4(x[n], E0, Einf, C, h) + noise;
  }

}

model{
  E0 ~ gamma(1/0.1, 1/(0.1*hc_med));
  e1 ~ beta(1, 1);
  log10_C ~ uniform(lb_log10_C, ub_log10_C);
  h ~ gamma(1.5, 0.5);
  s_lc ~ normal(0, 0.1);
  s_y ~ normal(0, 0.1);
  k ~ gamma(2, 0.2);

  for (n in 1:N) {
    pi[n] ~ beta(1, 19);
    target += log_mix(pi[n],
                      neg_binomial_2_lpmf(y[n] | u[n], 1/(5 + k)/s_y),
                      neg_binomial_2_lpmf(y[n] | u[n], 1/s_y));
  }

  for (n in 1:N_hc) {
    pi_hc[n] ~ beta(1, 19);
    target += log_mix(pi_hc[n],
                      neg_binomial_2_lpmf(hc[n] | E0 + noise, 1/(5 + k)/s_y),
                      neg_binomial_2_lpmf(hc[n] | E0 + noise, 1/s_y));
  }

  for (n in 1:N_lc) {
    lc[n] ~ neg_binomial_2(noise, 1/s_lc);
  }
}

generated quantities{
  // Expected target cell count
  real<lower=0> u_max;
  real<lower=0> u_min;
  array[N_new] real<lower=0> u_pred;
  // Predicted target cell count
  array[N_new] int<lower=0> y_pred;
  // Pointwise log likelihood
  vector[N + N_hc + N_lc] log_lik;

  u_max = E0 + noise;
  u_min = Einf + noise;

  // Predicted values
  for (n in 1:N_new) {
    u_pred[n] = logistic4(x_new[n], E0, Einf, C, h) + noise;
    y_pred[n] = neg_binomial_2_rng(u_pred[n], 1/s_y);
  }

  // Pointwise log likelihood
  for (n in 1:N) {
    log_lik[n] = log_mix(pi[n],
                         neg_binomial_2_lpmf(y[n] | u[n], 1/(5 + k)/s_y),
                         neg_binomial_2_lpmf(y[n] | u[n], 1/s_y));
  }
  for (n in 1:N_hc) {
    log_lik[N + n] = log_mix(pi_hc[n],
                             neg_binomial_2_lpmf(hc[n] | E0 + noise, 1/(5 + k)/s_y),
                             neg_binomial_2_lpmf(hc[n] | E0 + noise, 1/s_y));
  }
  for (n in 1:N_lc) {
    log_lik[N + N_hc + n] = neg_binomial_2_lpmf(lc[n] | noise, 1/s_lc);
  }
}

