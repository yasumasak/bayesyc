
functions{
  real logistic4(real d, real E0, real E1, real C, real h){
    return( E1 + (E0 - E1) / (1 + pow(d/C, h)) );
  }
}

data{
  // Observed variables
  int<lower=0> N; // The number of data points
  array[N] real<lower=0> y; // target cell viability
  array[N] real<lower=0> x; // Concentrations
  real<lower=0> E0; // Target cell viability without treatment

  // hyperparameters
  real lb_log10_C; // Lower bound for log10(C)
  real ub_log10_C; // Upper bound for log10(C)

  // Variables for prediction
  int<lower=0> N_new; // The number of data points
  array[N_new] real<lower=0> x_new; // Concentrations
}

parameters{
  // For dose response curve
  real<lower=0, upper=1> e1;
  real log10_C;
  real<lower=0> h;
  real<lower=0> s_y;

  // For outlier
  array[N] real<lower=0, upper=1> pi; // normal: 0, outlier: 1
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
    u[n] = logistic4(x[n], E0, Einf, C, h);
  }

}

model{
  e1 ~ beta(1, 1);
  log10_C ~ uniform(lb_log10_C, ub_log10_C);
  h ~ gamma(1.5, 0.5);
  s_y ~ normal(0, 0.1);
  k ~ gamma(2, 0.2);

  for (n in 1:N) {
    pi[n] ~ beta(1, 19);
    target += log_mix(pi[n],
                      gamma_lpdf(y[n] | 1/(5 + k)/s_y, 1/u[n]/(5 + k)/s_y),
                      gamma_lpdf(y[n] | 1/s_y, 1/u[n]/s_y));
  }

}

generated quantities{
  // Expected target cell viability
  real<lower=0> u_max;
  real<lower=0> u_min;
  array[N_new] real<lower=0> u_pred;
  // Predicted target cell viability
  array[N_new] real<lower=0> y_pred;
  // Pointwise log likelihood
  vector[N] log_lik;

  u_max = E0;
  u_min = Einf;

  // Predicted values
  for (n in 1:N_new) {
    u_pred[n] = logistic4(x_new[n], E0, Einf, C, h);
    y_pred[n] = gamma_rng(1/s_y, 1/u_pred[n]/s_y);
  }

  // Pointwise log likelihood
  for (n in 1:N) {
    log_lik[n] = log_mix(pi[n],
                         gamma_lpdf(y[n] | 1/(5 + k)/s_y, 1/u[n]/(5 + k)/s_y),
                         gamma_lpdf(y[n] | 1/s_y, 1/u[n]/s_y));
  }
}

