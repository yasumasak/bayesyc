
data{

  int<lower=0> N1; // The number of data points for drug1 only
  int<lower=0> N2; // The number of data points for drug2 only
  array[N1] int<lower=0> y1; // target cell count with drug1
  array[N2] int<lower=0> y2; // target cell count with drug2
  array[N1] real<lower=0> x1; // Concentrations for drug1
  array[N2] real<lower=0> x2; // Concentrations for drug2

  int<lower=0> Nh; // The number of high control data points
  array[Nh] int<lower=0> hc; // High control target cell count
  real<lower=0> hc_med; // Median of high control target cell count

  // hyperparameters
  real lb_log10_C; // Lower bound for log10(C)
  real ub_log10_C; // Upper bound for log10(C)

}


parameters{

  // For dose response curve
  real<lower=0> E0;
  real<lower=0, upper=1> e1;
  real<lower=0, upper=1> e2;
  real log10_C1;
  real log10_C2;
  real<lower=0> h1;
  real<lower=0> h2;
  real<lower=0> s_y;

  // For outlier
  array[N1] real<lower=0, upper=1> pi1; // normal: 0, outlier: 1
  array[N2] real<lower=0, upper=1> pi2; // normal: 0, outlier: 1
  array[Nh] real<lower=0, upper=1> pih; // normal: 0, outlier: 1
  real<lower=0> k;

}


transformed parameters{
  real<lower=0> E1;
  real<lower=0> E2;
  real<lower=0> C1;
  real<lower=0> C2;

  E1 = E0 * e1;
  E2 = E0 * e2;
  C1 = pow(10, log10_C1);
  C2 = pow(10, log10_C2);
}


model{

  // Expected target cell count
  array[N1] real u1;
  array[N2] real u2;

  log10_C1 ~ uniform(lb_log10_C, ub_log10_C);
  log10_C2 ~ uniform(lb_log10_C, ub_log10_C);
  h1 ~ gamma(1.5, 0.5);
  h2 ~ gamma(1.5, 0.5);
  E0 ~ gamma(1/0.1, 1/(0.1*hc_med));
  e1 ~ beta(1, 1);
  e2 ~ beta(1, 1);
  s_y ~ normal(0, 0.1);

  k ~ gamma(2, 0.2);

  for (n in 1:N1) {
    u1[n] = E1 + (E0 - E1) / (1 + pow(x1[n]/C1, h1));
    pi1[n] ~ beta(1, 19);
    target += log_mix(pi1[n],
                      neg_binomial_2_lpmf(y1[n] | u1[n], 1/(5 + k)/s_y),
                      neg_binomial_2_lpmf(y1[n] | u1[n], 1/s_y));
  }

  for (n in 1:N2) {
    u2[n] = E2 + (E0 - E2) / (1 + pow(x2[n]/C2, h2));
    pi2[n] ~ beta(1, 19);
    target += log_mix(pi2[n],
                      neg_binomial_2_lpmf(y2[n] | u2[n], 1/(5 + k)/s_y),
                      neg_binomial_2_lpmf(y2[n] | u2[n], 1/s_y));
  }

  for (n in 1:Nh) {
    pih[n] ~ beta(1, 19);
    target += log_mix(pih[n],
                      neg_binomial_2_lpmf(hc[n] | E0, 1/(5 + k)/s_y),
                      neg_binomial_2_lpmf(hc[n] | E0, 1/s_y));
  }

}


generated quantities{

  // Expected target cell count
  real<lower=0> u_max;
  real<lower=0> u_min;

  u_max = E0;
  u_min = E1;

}

