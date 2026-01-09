
functions{

  real musyc(real d1, real d2, real E0, real E1, real E2, real E3, real C1, real C2, real h1, real h2, real a12, real a21, real g12, real g21, real r1r, real r2r){
    real d1h1;
    real d2h2;
    real C1h1;
    real C2h2;
    real U;
    real A1;
    real A2;
    real r1;
    real r2;

    C1h1 = pow(C1,h1);
    C2h2 = pow(C2,h2);
    r1 = r1r/C1h1;
    r2 = r2r/C2h2;
    d1h1 = pow(d1,h1);
    d2h2 = pow(d2,h2);
    U=(r1*r2*pow((r1*C1h1),g21)*C1h1*C2h2+r1*r2*pow((r2*C2h2),g12)*C1h1*C2h2+pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)*C1h1+pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)*C2h2)/(d1h1*r1*r2*pow((r1*C1h1),g21)*C2h2+d1h1*r1*r2*pow((r2*C2h2),g12)*C2h2+d1h1*r1*pow(r2,(g12+1))*pow(a12*d2, g12*h2)*C2h2+d1h1*r1*pow(r2,g12)*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)+d1h1*pow(r1,(g21+1))*pow(r2,g12)*pow(a21*d1, g21*h1)*pow(a12*d2, g12*h2)+d1h1*pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)+d2h2*r1*r2*pow((r1*C1h1),g21)*C1h1+d2h2*r1*r2*pow((r2*C2h2),g12)*C1h1+d2h2*pow(r1,(g21+1))*r2*pow(a21*d1, g21*h1)*C1h1+d2h2*pow(r1,g21)*r2*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)+d2h2*pow(r1,g21)*pow(r2,(g12+1))*pow(a21*d1, g21*h1)*pow(a12*d2, g12*h2)+d2h2*pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)+r1*r2*pow((r1*C1h1),g21)*C1h1*C2h2+r1*r2*pow((r2*C2h2),g12)*C1h1*C2h2+pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)*C1h1+pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)*C2h2);
    A1=(d1h1*r1*r2*pow((r1*C1h1),g21)*C2h2+d1h1*r1*r2*pow((r2*C2h2),g12)*C2h2+d1h1*pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)+d2h2*pow(r1,g21)*r2*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12))/(d1h1*r1*r2*pow((r1*C1h1),g21)*C2h2+d1h1*r1*r2*pow((r2*C2h2),g12)*C2h2+d1h1*r1*pow(r2,(g12+1))*pow(a12*d2, g12*h2)*C2h2+d1h1*r1*pow(r2,g12)*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)+d1h1*pow(r1,(g21+1))*pow(r2,g12)*pow(a21*d1, g21*h1)*pow(a12*d2, g12*h2)+d1h1*pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)+d2h2*r1*r2*pow((r1*C1h1),g21)*C1h1+d2h2*r1*r2*pow((r2*C2h2),g12)*C1h1+d2h2*pow(r1,(g21+1))*r2*pow(a21*d1, g21*h1)*C1h1+d2h2*pow(r1,g21)*r2*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)+d2h2*pow(r1,g21)*pow(r2,(g12+1))*pow(a21*d1, g21*h1)*pow(a12*d2, g12*h2)+d2h2*pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)+r1*r2*pow((r1*C1h1),g21)*C1h1*C2h2+r1*r2*pow((r2*C2h2),g12)*C1h1*C2h2+pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)*C1h1+pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)*C2h2);
    A2=(d1h1*r1*pow(r2,g12)*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)+d2h2*r1*r2*pow((r1*C1h1),g21)*C1h1+d2h2*r1*r2*pow((r2*C2h2),g12)*C1h1+d2h2*pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21))/(d1h1*r1*r2*pow((r1*C1h1),g21)*C2h2+d1h1*r1*r2*pow((r2*C2h2),g12)*C2h2+d1h1*r1*pow(r2,(g12+1))*pow(a12*d2, g12*h2)*C2h2+d1h1*r1*pow(r2,g12)*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)+d1h1*pow(r1,(g21+1))*pow(r2,g12)*pow(a21*d1, g21*h1)*pow(a12*d2, g12*h2)+d1h1*pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)+d2h2*r1*r2*pow((r1*C1h1),g21)*C1h1+d2h2*r1*r2*pow((r2*C2h2),g12)*C1h1+d2h2*pow(r1,(g21+1))*r2*pow(a21*d1, g21*h1)*C1h1+d2h2*pow(r1,g21)*r2*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)+d2h2*pow(r1,g21)*pow(r2,(g12+1))*pow(a21*d1, g21*h1)*pow(a12*d2, g12*h2)+d2h2*pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)+r1*r2*pow((r1*C1h1),g21)*C1h1*C2h2+r1*r2*pow((r2*C2h2),g12)*C1h1*C2h2+pow(r1,(g21+1))*pow(a21*d1, g21*h1)*pow((r2*C2h2),g12)*C1h1+pow(r2,(g12+1))*pow(a12*d2, g12*h2)*pow((r1*C1h1),g21)*C2h2);

    return( U*E0 + A1*E1 + A2*E2 + (1-(U+A1+A2))*E3 );
  }

}


data{
  // Observed variables
  int<lower=0> Nc; // The number of data points for combination drug
  array[Nc] int yc; // target cell count with combination drug
  array[Nc] real<lower=0> dc1; // Concentrations for drug1 in combination
  array[Nc] real<lower=0> dc2; // Concentrations for drug2 in combination
  
  real<lower=0> hc_med; // Median of high control target cell count

  // Fixed parameters
  real<lower=0> r1r;
  real<lower=0> r2r;
  real<lower=0> gamma12;
  real<lower=0> gamma21;

  // Variables for prediction
  int<lower=0> Nc_new; // The number of data points
  array[Nc_new] real<lower=0> dc1_new; // Concentrations
  array[Nc_new] real<lower=0> dc2_new; // Concentrations
  
  // hyperparameters
  real lb_log10_C; // Lower bound for log10(C)
  real ub_log10_C; // Upper bound for log10(C)
}


parameters{

  // For basic dose response curve
  real<lower=0> E0;
  real<lower=0, upper=1> e1;
  real<lower=0, upper=1> e2;
  real<lower=0, upper=1> e3;
  real log10_C1;
  real log10_C2;
  real<lower=0> h1;
  real<lower=0> h2;
  real<lower=0> s_y;

  // For outlier
  array[Nc] real<lower=0, upper=1> pic; // normal: 0, outlier: 1
  real<lower=0> k;

  // For potency
  real log_alpha12;
  real log_alpha21;

}


transformed parameters{

  real<lower=0> C1;
  real<lower=0> C2;
  real<lower=0> E1;
  real<lower=0> E2;
  real<lower=0> E3;

  // Potency
  real<lower=0> alpha12;
  real<lower=0> alpha21;

  // Efficacy
  real beta1;
  real beta2;

  // Expected target cell count
  array[Nc] real uc;

  C1 = pow(10, log10_C1);
  C2 = pow(10, log10_C2);

  E1 = E0 * e1;
  E2 = E0 * e2;
  E3 = E0 * e3;

  alpha12 = exp(log_alpha12);
  alpha21 = exp(log_alpha21);

  beta1 = (E1 - E3)/(E0 - E1);
  beta2 = (E2 - E3)/(E0 - E2);

  for (n in 1:Nc) {
    uc[n] = musyc(dc1[n], dc2[n], E0, E1, E2, E3, C1, C2, h1, h2, alpha12, alpha21, gamma12, gamma21, r1r, r2r);
  }
}


model{

  E0 ~ gamma(1/0.1, 1/(0.1*hc_med));
  log10_C1 ~ uniform(lb_log10_C, ub_log10_C);
  log10_C2 ~ uniform(lb_log10_C, ub_log10_C);
  h1 ~ gamma(1.5, 0.5);
  h2 ~ gamma(1.5, 0.5);
  e1 ~ beta(1, 1);
  e2 ~ beta(1, 1);
  e3 ~ beta(1, 1);
  s_y ~ normal(0, 0.1); // CV ^ 2

  log_alpha12 ~ normal(0, 5);
  log_alpha21 ~ normal(0, 5);

  k ~ gamma(2, 0.2);

  for (n in 1:Nc) {
    pic[n] ~ beta(1, 19);
    target += log_mix(pic[n],
                      neg_binomial_2_lpmf(yc[n] | uc[n], 1/(5 + k)/s_y),
                      neg_binomial_2_lpmf(yc[n] | uc[n], 1/s_y));
  }
}


generated quantities{

  // Expected target cell count
  array[Nc_new] real<lower=0> uc_pred;
  real<lower=0> u_max;

  // Predicted target cell count
  array[Nc_new] int yc_pred;

  // Pointwise log likelihood
  vector[Nc] log_lik;

  u_max = E0;

  for (n in 1:Nc_new) {
    uc_pred[n] = musyc(dc1_new[n], dc2_new[n], E0, E1, E2, E3, C1, C2, h1, h2, alpha12, alpha21, gamma12, gamma21, r1r, r2r);
    yc_pred[n] = neg_binomial_2_rng(uc_pred[n], 1/s_y);
  }

  for (n in 1:Nc) {
    log_lik[n] = log_mix(pic[n],
                         neg_binomial_2_lpmf(yc[n] | uc[n], 1/(5 + k)/s_y),
                         neg_binomial_2_lpmf(yc[n] | uc[n], 1/s_y));
  }

}

