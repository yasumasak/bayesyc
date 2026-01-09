
#' Get Gamma distribution parameters from MCMC samples
#'
#' @param mx_samples Matrix of MCMC samples
#' @param rv random variable of interest 
#' @param coef_var coefficient applied to variance
#'
#' @return A vector of gamma distribution parameters, alpha and beta
get_gamma_param = function(mx_samples, rv, coef_var=1){
  if(any(colnames(mx_samples) == rv)){
    samples <- mx_samples[, rv]
    mean <- mean(samples)
    var <- var(samples) * coef_var
    alpha <- (mean^2) / var
    beta <- mean / var
    return(c(alpha, beta))
  }else{
    return(NULL)
  }
}


#' Get Beta distribution parameters from MCMC samples
#'
#' @param mx_samples Matrix of MCMC samples
#' @param rv random variable of interest 
#' @param coef_var coefficient applied to variance
#'
#' @return A vector of beta distribution parameters, alpha and beta
get_beta_param = function(mx_samples, rv, coef_var=1){
  if(any(colnames(mx_samples) == rv)){
    samples <- mx_samples[, rv]
    mean <- mean(samples)
    var <- var(samples) * coef_var
    alpha <- mean^2 * (1 - mean) / var - mean
    beta <- alpha / mean - alpha
    return(c(alpha, beta))
  }else{
    return(NULL)
  }
}