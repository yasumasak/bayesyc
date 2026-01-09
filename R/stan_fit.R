
#' @title Basic Stan fit information class
#'
#' @description
#' This class behaves as an abstract class as \code{estimates}, \code{predictions} and \code{stats} are defined by subclass,
#' using 'Bridge pattern'.  
#' This class holds Stan fitting result object, input list and other basic information of Stan fit result.
#'
#' @name Fit_Info
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{fit}: }{fit (stanfit or CmdStanMCMC) object retured from MCMC sampling.}
#'   \item{\code{input}: }{list used as a input of MCMC sampling.}
#'   \item{\code{model}: }{model (stanmodel or CmdStanModel) object used for fitting.}
#'   \item{\code{par_dims}: }{dimensions of parameters in the model.}
#'   \item{\code{fit_summary}: }{matrix of summary stats of estimated parameters.}
#'   \item{\code{thre_rhat}: }{numeric threshold for Rhat, used for judging success of fit.}
#'   \item{\code{sf_tc}: }{numeric scaling factor used for T/C calculation.}
#'   \item{\code{thre_tc}: }{numeric threshold for T/C in %, used for judging IC50.}
#'   \item{\code{thre_ic50}: }{numeric threshold: display IC50 if within min_dose/\code{thre_ic50} to max_dose*\code{thre_ic50}.}
#'   \item{\code{success_fit}: }{TRUE if successful fit, otherwise FALSE.}
#'   \item{\code{ok_Rhat}: }{TRUE if Rhat of all parameters < \code{thre_rhat}, otherwise FALSE.}
#'   \item{\code{ok_TC}: }{TRUE if minimum T/C of predicted u < \code{thre_tc}, otherwise FALSE.}
#'   \item{\code{ok_IC50}: }{TRUE if minimum T/C of predicted u < \code{thre_tc}, otherwise FALSE.}
#'   \item{\code{waic}: }{matrix of estimates from loo::waic.}
#'   \item{\code{looic}: }{matrix of estimates from loo::loo.}
#'   \item{\code{estimates}: }{Estimates class object which is defined by subclass.}
#'   \item{\code{predictions}: }{Predictions class object which is defined by subclass.}
#'   \item{\code{stats}: }{Stats class object which is defined by subclass.}
#' }
#' @keywords fit, Stan, base
#' @family Base
Fit_Info <- R6::R6Class("Fit_Info",

  public = list(
    fit = NULL, 
    input = NULL, 
    model = NULL,
    par_dims = NULL,
    fit_summary = NULL,
    thre_rhat = NULL, sf_tc = NULL, thre_tc = NULL, thre_ic50 = NULL,
    success_fit = NULL, ok_Rhat = NULL, ok_TC = NULL, ok_IC50 = NULL,
    waic = NULL, looic = NULL,
    # estimates, stats, predictions are set by subclass
    estimates = NULL,
    stats = NULL,
    predictions = NULL,

    #' @description
    #' Set \code{fit} and \code{input}, \code{fit_summary}, and \code{sf_tc} 
    #' using information from \code{fit}.
    #' \code{sf_tc} is set to 100 / median of estimated u_max as a default.
    #' @param fit fitting result object.
    #' @param input list of input data used for Stan sampling.
    #' @param model model used for fitting.
    initialize = function(fit, input, model, ...){
      self$input <- input
      self$model <- model
      private$set_fit(fit)
      
      if(!is.null(self$fit)){
        if(class(self$model)[1] == "CmdStanModel"){
          tbl_summary <- self$fit$summary(variables=NULL, 
                                          quantile, .args=list(probs=c(0.025, .5, .975)), 
                                          "mean", Rhat="rhat", "ess_bulk", "ess_tail")
          self$fit_summary <- as.matrix(dplyr::select(tbl_summary, -variable))
          rownames(self$fit_summary) <- tbl_summary$variable
          self$par_dims <- self$fit$metadata()$stan_variable_sizes
        }else{
          self$fit_summary <- rstan::summary(self$fit)$summary
          self$par_dims <- self$fit@par_dims
        }
        self$sf_tc <- 100 / self$fit_summary["u_max", "50%"]
      }else{
        self$sf_tc <- 100 / self$input$hc_med
      }
    },

    #' @description
    #' Set \code{sf_tc} (scaling factor of T/C).
    set_sf_tc = function(sf_tc){
      self$sf_tc <- sf_tc
    },
    
    #' @description
    #' Check success of fit.
    #' @param thre_rhat threshold for Rhat
    check_fit = function(thre_rhat){
      self$thre_rhat <- thre_rhat
      if(!is.null(self$fit)){
        # Check success of fit
        if(!all(self$fit_summary[,"Rhat"] < self$thre_rhat, na.rm=T)){
          self$success_fit <- FALSE
          self$ok_Rhat <- FALSE
        }else{
          self$success_fit <- TRUE
          self$ok_Rhat <- TRUE
        }
      }else{
        self$success_fit <- FALSE
        self$ok_Rhat <- NA
      }
    },
    
    #' @description
    #' Set \code{thre_tc} and \code{$ok_TC}.
    #' @param thre_tc threshold for minimum T/C.
    check_tc = function(thre_tc){
      self$thre_tc <- thre_tc
      if(!is.null(self$fit)){
        if(!is.null(self$predictions)){
          self$ok_TC <- self$predictions$check_tc(self$thre_tc)
        }else{
          warning("T/C can not be checked. Do 'set_predictions()' beforehand.")
        }
      }else{
        self$ok_TC <- NA
      }
    },

    #' @description
    #' Set \code{thre_tc} and \code{$ok_TC}.
    #' @param thre_tc threshold for detected IC50.
    check_ic50 = function(thre_ic50){
      self$thre_ic50 <- thre_ic50
      if(!is.null(self$fit)){
        if(!is.null(self$ok_TC)){
          self$ok_IC50 <- self$ok_TC & self$estimates$check_ic50(self$input, self$thre_ic50)
        }else{
          warning("IC50 can not be checked. Do 'check_tc()' beforehand.")
        }
      }else{
        self$ok_IC50 <- NA
      }
    },
    
    #' @description
    #' Return pi[] variable names if estimated.
    #' @return character "pi[]"
    get_pi_names = function(){
      if(!is.null(self$par_dims)){
        if(!is.null(self$par_dims$pi)){
          return(paste0("pi[", 1:(self$par_dims$pi),"]"))
        }
      }
    },

    #' @description
    #' Calculate and set \code{waic}, \code{looic} and return Estimate ± SE of them.
    calculate_IC = function(){
      if(!is.null(self$fit)){
        # ! This function needs vector log_lik variable in Stan code.
        if(class(self$model)[1] == "CmdStanModel"){
          mx_log_lik <- self$fit$draws(variables="log_lik", format="draws_matrix")
        }else{
          mx_log_lik <- loo::extract_log_lik(self$fit)
        }
        # WAIC
        self$waic <- loo::waic(mx_log_lik)$estimates
        waic_str <- paste0("WAIC: ", self$waic["waic", "Estimate"], " ± ", self$waic["waic", "SE"])
        # LOO
        self$looic <- loo::loo(mx_log_lik)$estimates
        looic_str <- paste0("LOO: ", self$looic["looic", "Estimate"], " ± ", self$looic["looic", "SE"])
        return(list(WAIC=waic_str, LOO=looic_str))
      }
    },
    
    #' @description
    #' Return c(\code{success_fit}, \code{ok_Rhat}, \code{ok_TC}, \code{ok_IC50}) for stats.
    get_judgement = function(){
      v_judgement <- c(self$success_fit, self$ok_Rhat, self$ok_TC)
      names(v_judgement) <- c("Success of fit", paste0("Rhat < ",self$thre_rhat), 
                              paste0("min(T/C) < ",self$thre_tc,"%"))
      if(!is.null(self$ok_IC50)){
        v_judgement <- append(v_judgement, c(`IC50 detected` = self$ok_IC50))
      }
      return(v_judgement)
    },
    
    #' @description
    #' Choose one Fit_Info object from \code{self} and \code{fit_info_alt} by comparing their WAIC.
    #' Myself (\code{self}) is the primary choice if alt_first=FALSE (default).
    #' Alternative (\code{fit_info_alt}) is selected if it satisfies the criteria.
    #' @param fit_info_alt an alternative Fit_Info object.
    #' @param thre_rhat threshold for Rhat
    #' @param alt_first FLASE (default). TRUE to prioritize alternative Fit_Info object.
    #' @return A selected Fit_Info object.
    select_fit_info = function(fit_info_alt, thre_rhat, alt_first=FALSE){
      # Check fit
      self$check_fit(thre_rhat)
      fit_info_alt$check_fit(thre_rhat)
      # Calculate IC
      self$calculate_IC()
      fit_info_alt$calculate_IC()
      
      # Check success of fit
      if(!fit_info_alt$success_fit){
        return(self)
      }else if(!self$success_fit){
        return(fit_info_alt)
      }
      
      # Check WAIC when both models fit successfully
      if(alt_first){
        # Prioritize alternative Fit_Info object
        condition <- (fit_info_alt$waic["waic", "Estimate"] < self$waic["waic", "Estimate"] + self$waic["waic", "SE"])
      }else{
        # Prioritize self Fit_Info object
        condition <- (fit_info_alt$waic["waic", "Estimate"] + fit_info_alt$waic["waic", "SE"] < self$waic["waic", "Estimate"])
      }
      if(condition){
        return(fit_info_alt)
      }else{
        return(self)
      }
    },
    
    #' @description
    #' Return a matrix of MCMC draws.
    #' @return matrix of MCMC draws
    get_draws = function(){
      if(class(self$model)[1] == "CmdStanModel"){
        draws <- self$fit$draws(format = "draws_matrix") 
      }else{
        df_draws <- as.data.frame(rstan::extract(self$fit))
        colnames(df_draws) <- gsub("\\.(\\d+)$", "[\\1]", colnames(df_draws))
        draws <- as.matrix(df_draws)
      }
      return(draws)
    }
  ),
  
  private = list(
    #' @description
    #' Set \code{fit} if Stan fitting result object is valid.
    #' @param fit fitting result object
    set_fit = function(fit){
      if(!is.null(fit)){
        if(class(self$model)[1] == "stanmodel"){
          if(length(fit@sim) > 0) self$fit <- fit
        }else{
          self$fit <- fit
        }
      }
    }
  )
)

#' @title Mono-phasic Logistic4 fit information class
#'
#' @description
#' Concrete class of Fit_Info for  Mono-phasic Logistic4 fit.
#'
#' @name Fit_MonoPhasic4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Fit_Info}}{}
#' }
#' 
#' @section Public fields:
#' 
#' @keywords fit, Stan
#' @family Logistic4
#' 
#' @export
Fit_MonoPhasic4 <- R6::R6Class("Fit_MonoPhasic4", inherit = Fit_Info,
                               
  public = list(
    
    #' @description
    #' Run \code{initialize} declared in the superclass.
    #' Set instance of \code{Estimates_MonoPhasic4} to \code{estimates} declared in the superclass.
    #' Set instance of \code{Stats_Logistic4} to \code{stats} declared in the superclass.
    #' @param fit fitting result object.
    #' @param input list used as a input of Stan sampling.
    #' @param model model used for fitting.
    initialize = function(fit, input, model){
      super$initialize(fit, input, model)
      self$estimates <- Estimates_MonoPhasic4$new(self)
      self$stats <- Stats_Logistic4$new(self)
    },
    
    #' @description
    #' Set instance of \code{Predictions_Curve} to \code{predictions} declared in the superclass.
    set_predictions = function(){
      self$predictions <- Predictions_Curve$new(self)
    }
  )
)

#' @title Bi-phasic Logistic4 fit information class
#'
#' @description
#' Concrete class of Fit_Info for Bi-phasic Logistic4 fit.
#' Biphasic model: Di Veroli et al. Scientific Reports volume 5, Article number: 14701, 2015. 
#' https://www.nature.com/articles/srep14701
#'
#' @name Fit_BiPhasic4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Fit_Info}}{}
#' }
#' 
#' @section Public fields:
#' 
#' @keywords fit, Stan
#' @family Logistic4
#' 
#' @export
Fit_BiPhasic4 <- R6::R6Class("Fit_BiPhasic4", inherit = Fit_Info,
                             
  public = list(

    #' @description
    #' Run \code{initialize} declared in the superclass.
    #' Set instance of \code{Estimates_BiPhasic4} to \code{estimates} declared in the superclass.
    #' Set instance of \code{Stats_Logistic4} to \code{stats} declared in the superclass.    
    #' @param fit fitting result object.
    #' @param input list used as a input of Stan sampling.
    #' @param model model used for fitting.
    initialize = function(fit, input, model){
      super$initialize(fit, input, model)
      self$estimates <- Estimates_BiPhasic4$new(self)
      self$stats <- Stats_Logistic4$new(self)
    },

    #' @description
    #' Set instance of \code{Predictions_Curve} to \code{predictions} declared in the superclass.    
    set_predictions = function(){
      self$predictions <- Predictions_Curve$new(self)
    }
  )
)

#' @title Fit_Dual fit information class
#'
#' @description
#' Concrete class of Fit_Info for Fit_Dual fit.
#'
#' @name Fit_Dual
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Fit_Info}}{}
#' }
#' 
#' @section Public fields:
#' 
#' @keywords fit, Stan
#' @family Logistic4
#' 
#' @export
Fit_Dual <- R6::R6Class("Fit_Dual", inherit = Fit_Info,
                          
  public = list(

    #' @description
    #' Run \code{initialize} declared in the superclass.
    #' @param fit fitting result object.
    #' @param input list used as a input of Stan sampling.
    #' @param model model used for fitting.
    initialize = function(fit, input, model){
      super$initialize(fit, input, model)
    },

    #' @description
    #' Return estimated parameters from Posterior distributions of Dual for BayeSyC
    get_bayesyc_params = function(){
      # Parameters from Posterior distributions of Dual
      draws_d <- super$get_draws()
      E0_param <- get_gamma_param(draws_d, "E0")
      e1_param <- get_beta_param(draws_d, "e1")
      e2_param <- get_beta_param(draws_d, "e2")
      C1_param <- get_gamma_param(draws_d, "C1")
      C2_param <- get_gamma_param(draws_d, "C2")
      h1_param <- get_gamma_param(draws_d, "h1")
      h2_param <- get_gamma_param(draws_d, "h2")
      s_y_mean <- mean(draws_d[,"s_y"])
      
      # Input data for BayeSyC
      return(list(E0_alpha=E0_param[1], E0_beta=E0_param[2],
                  e1_alpha=e1_param[1], e1_beta=e1_param[2],
                  e2_alpha=e2_param[1], e2_beta=e2_param[2],
                  C1_alpha=C1_param[1], C1_beta=C1_param[2],
                  C2_alpha=C2_param[1], C2_beta=C2_param[2],
                  h1_alpha=h1_param[1], h1_beta=h1_param[2],
                  h2_alpha=h2_param[1], h2_beta=h2_param[2],
                  s_y=s_y_mean))
    }
  )
)

#' @title Fit_BayeSyC4 fit information class
#'
#' @description
#' Concrete class of Fit_Info for Fit_BayeSyC4 fit.
#'
#' @name Fit_BayeSyC4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Fit_Info}}{}
#' }
#' 
#' @section Public fields:
#' 
#' @keywords fit, Stan
#' @family BayeSyC4
#' 
#' @export
Fit_BayeSyC4 <- R6::R6Class("Fit_BayeSyC4", inherit = Fit_Info,
                          
  public = list(

    #' @description
    #' Run \code{initialize} declared in the superclass.
    #' Set instance of \code{Estimates_BayeSyC4} to \code{estimates} declared in the superclass.
    #' Set instance of \code{Stats_BayeSyC4} to \code{stats} declared in the superclass.
    #' @param fit fitting result object.
    #' @param input list used as a input of Stan sampling.
    #' @param model model used for fitting.
    initialize = function(fit, input, model){
      super$initialize(fit, input, model)
      self$estimates <- Estimates_BayeSyC4$new(self)
      self$stats <- Stats_BayeSyC4$new(self)
    },

    #' @description
    #' Set instance of \code{Predictions_Surface} to \code{predictions} declared in the superclass.      
    set_predictions = function(){
      self$predictions <- Predictions_Surface$new(self)
    }
  )
)


#' @title Basic Estimates class
#'
#' @description
#' This class behaves as an abstract class as \code{est_params} is defined by subclass.
#' This class holds \code{df} that contains estimated values.
#'
#' @name Estimates
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{est_params}: }{parameter names to set \code{df}.}
#'   \item{\code{df}: }{data.frame containing estimated values.
#'   Rows are \code{est_params} and columns are c("2.5\%", "50\%", "97.5\%").}
#' }
#' @keywords estimates, Stan, base
#' @family Base
Estimates <- R6::R6Class("Estimates", 
                         
  public = list(
    est_params = NULL,
    df = NULL,

    #' @description
    #' Set \code{df} containing estimated values.
    #' Assumed to be called from subclass with the argument \code{fit_summary}.
    #' @param fit_summary matrix of summary stats of estimated parameters.
    set_estimates = function(fit_summary){
      if(!is.null(fit_summary)){
        self$est_params <- self$get_valid_params(fit_summary, self$est_params)
        if(length(self$est_params) > 1){
          df_estimate <- transform(data.frame(signif(fit_summary[self$est_params, c("2.5%", "50%", "97.5%")], 3), check.names=F),
                                   parameter=self$est_params, check.names=F)
        }else if(length(self$est_params) == 1){
          df_estimate <- transform(data.frame(t(signif(fit_summary[self$est_params, c("2.5%", "50%", "97.5%")], 3)), check.names=F),
                                   parameter=self$est_params, check.names=F)
        }else{
          stop("Fetal error in 'set_estimates' function.")
        }
        rownames(df_estimate) <- df_estimate$parameter
        self$df <- df_estimate
      }
    },
    
    #' @description
    #' This method check the existence of \code{params_to_check} in \code{fit_summary},
    #' and return valid params
    #' @param fit_summary matrix of summary stats of estimated parameters.
    #' @param params_to_check parameters to be checked.
    get_valid_params = function(fit_summary, params_to_check){
      # Check if est_params exist in fit_summary
      valid_params <- NULL
      valid_names <- NULL
      for(i in seq(params_to_check)){
        if(any(rownames(fit_summary) == params_to_check[i])){
          valid_params <- c(valid_params, params_to_check[i])
          valid_names <- c(valid_names, names(params_to_check)[i])
        }
      }
      names(valid_params) <- valid_names
      return(valid_params)
    },
    
    #' @description
    #' Check IC50
    #' @param input input list object used for MCMC sampling.
    #' @param thre_ic50 numeric threshold for concentration range.
    check_ic50 = function(input, thre_ic50){
      # Check median IC50
      median_ic50 <- self$df["C", "50%"]
      if(!is.na(median_ic50)){
        if((min(input$x) / thre_ic50 < median_ic50) & (median_ic50 < thre_ic50 * max(input$x))){
          return(TRUE)
        }else{
          return(FALSE)
        }
      }else{
        return(NA)
      }
    }
  )
)

#' @title Estimates class for monophasic logistic curve model
#'
#' @description
#' Concrete class of Estimates for monophasic logistic curve model.
#'
#' @name Estimates_MonoPhasic4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Estimates}}{}
#' }
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{phase}: }{number of phases in logistic curve. Set "mono" in this class.}
#'   \item{\code{representative}: }{representative phase in multiphasic model. Set NA in this class.}
#' }
#' @keywords estimates, Stan
#' @family Logistic4
Estimates_MonoPhasic4 <- R6::R6Class("Estimates_MonoPhasic4",  inherit = Estimates,
                                     
  public = list(
    phase = "mono",
    representative = NA,

    #' @description
    #' Set \code{est_params} and \code{df} declared in the superclass 
    #' using the argument \code{fit_info}.
    #' @param fit_info Fit_info class object.
    initialize = function(fit_info){
      self$est_params <- c("u_max", "u_min", "noise", "E0", "Einf", "h", "C", "s_y", "s_lc")
      self$est_params <- c(self$est_params, fit_info$get_pi_names())
      super$set_estimates(fit_info$fit_summary)
    }
  )
)

#' @title Estimates class for biphasic logistic curve model
#'
#' @description
#' Concrete class of Estimates for biphasic logistic curve model.
#' Biphasic model: Di Veroli et al. Scientific Reports volume 5, Article number: 14701, 2015. 
#' https://www.nature.com/articles/srep14701
#'
#' @name Estimates_BiPhasic4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Estimates}}{}
#' }
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{phase}: }{number of phases in logistic curve. Set "bi" in this class.}
#'   \item{\code{representative}: }{Representative phase in multiphasic model.}
#' }
#' @keywords estimates, Stan
#' @family Logistic4
Estimates_BiPhasic4 <- R6::R6Class("Estimates_BiPhasic4",  inherit = Estimates,
                                   
  public = list(
    phase = "bi",
    representative = NULL,

    #' @description
    #' Set \code{est_params} and \code{df} declared in the superclass 
    #' using the argument \code{fit_info}.
    #' @param fit_info Fit_info class object.
    initialize = function(fit_info){
      self$est_params <- c("u_max", "u_min", "noise", "E0", "E1", "E2", "h1", "h2", "C1", "C2", "s_y", "s_lc")
      self$est_params <- c(self$est_params, fit_info$get_pi_names())
      super$set_estimates(fit_info$fit_summary)
      private$set_representative()
    }
  ),

  private = list(
    #' @description
    #' Set \code{self$representative}
    set_representative = function(){
      if(!is.null(self$df)){
        E0 <- self$df["u_max", "50%"]
        E1 <- self$df["E1", "50%"]
        E2 <- self$df["E2", "50%"]
        if((E0 - E1) > (E1 - E1 * E2 / E0)){
          self$df["h", ] <- self$df["h1", ]
          self$df["C", ] <- self$df["C1", ]
          self$representative <- "1"
        }else{
          self$df["h", ] <- self$df["h2", ]
          self$df["C", ] <- self$df["C2", ]
          self$representative <- "2"
        }
      }else{
        self$representative <- NA
      }
    }
  )
)

#' @title Estimates class for BayeSyC model
#'
#' @description
#' Concrete class of Estimates for BayeSyC model.
#'
#' @name Estimates_BayeSyC4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Estimates}}{}
#' }
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{beta_name}: }{"beta1" or "beta2" used to estimate efficacy.}
#' }
#' @keywords estimates, Stan
#' @family BayeSyC4
Estimates_BayeSyC4 <- R6::R6Class("Estimates_BayeSyC4",  inherit = Estimates,
                                
  public = list(
    beta_name = NULL,

    #' @description
    #' Set \code{est_params} and \code{df} declared in the superclass 
    #' using the argument \code{fit_info}.
    #' @param fit_info Fit_info class object.
    initialize = function(fit_info){
      self$est_params <- c("log_alpha12", "log_alpha21", "log_alpha", "alpha12", "alpha21", "alpha", "C1", "C2", "beta1", "beta2", "E0", "E1", "E2", "E3", "h1", "h2", "s_y", "noise", "s_lc")
      super$set_estimates(fit_info$fit_summary)
      private$set_beta() # Set beta_name "beta1" or "beta2"
    }
  ),

  private = list(
    #' @description
    #' Set \code{self$beta_name} "beta1" or "beta2".
    set_beta = function(){
      if(!is.null(self$df)){
        if(self$df["E1", "50%"] < self$df["E2", "50%"]){
          self$beta_name <- "beta1"
        }else{
          self$beta_name <- "beta2"
        }
        self$df["beta", ] <- self$df[self$beta_name, ]
      }else{
        self$beta_name <- "beta"
      }
    }
  )
)


#' @title Basic Stats class
#'
#' @description
#' This class behaves as an abstract class as \code{stats_params} is defined by subclass.
#' This class holds stats vector \code{values} that contains information of 
#' \code{stats_param} renamed as \code{stats_name}.
#'
#' @name Stats
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{stats_params}: }{parameter names to construct stats vector \code{values}.}
#'   \item{\code{values}: }{vector of stats that are 2.5\%, 50\%, 97.5\%, Rhat and Neff values of \code{stats_params}.}
#' }
#' @keywords estimates, Stan, base
#' @family Base
Stats <- R6::R6Class("Stats",
                     
  public = list(
    stats_params = NULL,
    values = NULL,

    #' @description
    #' Set stats vector \code{values}.
    #' Assumed to be called from subclass with the argument \code{fit_summary} and \code{estimates}.
    #' @param fit_summary matrix of summary stats of estimated parameters.
    #' @param estimates Estimates class object.
    #' @param model model used for fitting.
    set_stats = function(fit_summary, estimates, model){
      v_to_show <- c("2.5%", "50%", "97.5%")
      # Output summary for effective sample size is different between Stan and CmdStanR
      if(class(model)[1] == "CmdStanModel"){ v_eff <- c("ess_bulk", "ess_tail") }else{ v_eff <- "n_eff" }
      if(!is.null(fit_summary)){
        # With valid parameters 
        self$stats_params <- estimates$get_valid_params(estimates$df, self$stats_params)
        # Values to output
        v_stats <- as.character(t(signif(estimates$df[self$stats_params, v_to_show], 3)))
        v_param <- as.character(estimates$df[self$stats_params, "parameter"])
        v_n_eff <- as.character(round(fit_summary[v_param, v_eff]))
        v_rhat <- as.character(signif(fit_summary[v_param, "Rhat"], 3))
      }else{
        # NAs for no result
        v_stats <- rep(NA, length(self$stats_params) * length(v_to_show))
        v_n_eff <- rep(NA, length(self$stats_params) * length(v_eff))
        v_rhat <- rep(NA, length(self$stats_params))
      }
      # Names
      names(v_stats) <- paste(rep(names(self$stats_params), each=length(v_to_show)), v_to_show)
      names(v_rhat) <- paste("Rhat", names(self$stats_params))
      names(v_n_eff) <- paste(rep(v_eff, each=length(self$stats_params)), names(self$stats_params))
      self$values <- c(v_stats, v_rhat, v_n_eff)
    }
  )
)

#' @title Stats class for logistic curve model
#'
#' @description
#' Concrete class of Stats for logistic curve model.
#'
#' @name Stats_Logistic4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Stats}}{}
#' }
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{stats_model_info}: }{header of model info: c("Phase", "Rep-param").}
#' }
#' @keywords stats, Stan
#' @family Logistic4
Stats_Logistic4 <- R6::R6Class("Stats_Logistic4", inherit = Stats,
                               
  public = list(
    stats_model_info = c("Phase", "Rep-param"),

    #' @description
    #' Set \code{stats_params} and \code{values} declared in the superclass 
    #' using the argument \code{fit_info}.
    #' @param fit_info Fit_info class object.
    initialize = function(fit_info){
      self$stats_params <- c("C", "u_max", "u_min", "h", "s_y", "noise", "s_lc")
      names(self$stats_params) <- c("C", "E0", "Einf", "h", "s_y", "noise", "s_lc")
      super$set_stats(fit_info$fit_summary, fit_info$estimates, fit_info$model)
      # Logistic4 specific
      private$append_model_info(fit_info$estimates)
    }
  ),

  private = list(
    #' @description
    #' Append model info to \code{values} declared in the superclass.
    #' @param estimates Estimates class object.
    append_model_info = function(estimates){
      v_model <- c(estimates$phase, estimates$representative)
      names(v_model) <- self$stats_model_info
      self$values <- c(v_model, self$values)
    }
  )
)

#' @title Stats class for BayeSyC model
#'
#' @description
#' Concrete class of Stats for BayeSyC model.
#'
#' @name Stats_BayeSyC4
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Stats}}{}
#' }
#' 
#' @section Public fields:
#'
#' @keywords stats, Stan
#' @family BayeSyC4
Stats_BayeSyC4 <- R6::R6Class("Stats_BayeSyC4", inherit = Stats,
                            
  public = list(
                             
    #' @description
    #' Set \code{stats_params} and \code{values} declared in the superclass 
    #' using the argument \code{fit_info}.
    #' @param fit_info Fit_info class object.
    initialize = function(fit_info){
      self$stats_params <- c("alpha12", "alpha21", "alpha", "C1", "C2", fit_info$estimates$beta_name, "E0", "E1", "E2", "E3", "h1", "h2", "s_y", "noise", "s_lc")
      names(self$stats_params) <- c("Alpha12", "Alpha21", "Alpha", "C1", "C2", "Beta", "E0", "E1", "E2", "E3", "h1", "h2", "s_y", "noise", "s_lc")
      super$set_stats(fit_info$fit_summary, fit_info$estimates, fit_info$model)
    }
  )
)


#' @title Basic Predictions class
#'
#' @description
#' This class behaves as an abstract class as \code{y_pred_names} and \code{u_pred_names} 
#' are defined by subclass.
#' This class holds \code{df} that contains predicted values of \code{y_pred_names} and \code{u_pred_names}.
#'
#' @name Predictions
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{y_pred_names}: }{parameter names of predicted y to set \code{df}.}
#'   \item{\code{u_pred_names}: }{parameter names of predicted u to set \code{df}.}
#'   \item{\code{df}: }{data.frame containing predicted values that are 2.5, 50 and 97.5 percentile values 
#'   of \code{y_pred_names}, \code{u_pred_names} and their T/C values.}
#'   \item{\code{xx}: }{vector of new x values used in the prediction.}
#' }
#' @keywords predictions, Stan, base
#' @family Base
Predictions <- R6::R6Class("Predictions", 
                           
  public = list(
    y_pred_names = NULL,
    u_pred_names = NULL,
    df = NULL,
    xx = NULL,
    
    #' @description
    #' Set of predicted 2.5, 50 and 97.5 percentile values of y, u, tc_y and tc_u.
    #' Assumed to be called from subclass.
    #' @param fit_summary matrix of summary stats of estimated parameters.
    #' @param input list used as a input of Stan sampling.
    #' @param sf_tc numeric scaling factor to calculate T/C.
    set_predictions = function(fit_summary, input, sf_tc){
      if(!is.null(fit_summary)){
        # Predicted values of cells
        mx_y_pred <- fit_summary[self$y_pred_names, c("50%","2.5%","97.5%")]
        v_y_median <- signif(mx_y_pred[,"50%"], 3)
        v_y_lower <- signif(mx_y_pred[,"2.5%"], 3)
        v_y_upper <- signif(mx_y_pred[,"97.5%"], 3)
        # Expected values of cells
        mx_u_pred <- fit_summary[self$u_pred_names, c("50%","2.5%","97.5%")]
        v_u_median <- signif(mx_u_pred[,"50%"], 3)
        v_u_lower <- signif(mx_u_pred[,"2.5%"], 3)
        v_u_upper <- signif(mx_u_pred[,"97.5%"], 3)
        # T/C
        v_tc_y_median <- signif(mx_y_pred[,"50%"] * sf_tc, 3)
        v_tc_y_lower <- signif(mx_y_pred[,"2.5%"] * sf_tc, 3)
        v_tc_y_upper <- signif(mx_y_pred[,"97.5%"] * sf_tc, 3)
        v_tc_u_median <- signif(mx_u_pred[,"50%"] * sf_tc, 3)
        v_tc_u_lower <- signif(mx_u_pred[,"2.5%"] * sf_tc, 3)
        v_tc_u_upper <- signif(mx_u_pred[,"97.5%"] * sf_tc, 3)

        self$xx <- input$x_new
        self$df <- data.frame(y_median=v_y_median, y_lower=v_y_lower, y_upper=v_y_upper,
                     u_median=v_u_median, u_lower=v_u_lower, u_upper=v_u_upper,
                     tc_y_median=v_tc_y_median, tc_y_lower=v_tc_y_lower, tc_y_upper=v_tc_y_upper,
                     tc_u_median=v_tc_u_median, tc_u_lower=v_tc_u_lower, tc_u_upper=v_tc_u_upper)
      }
    },

    #' @description
    #' Return TRUE if any median T/Cs > \code{thre_tc} otherwise FALSE.
    #' @param thre_tc \% threshold for T/C 
    #' @return TRUE if minimum of T/C > \code{thre_tc} otherwise FALSE.
    check_tc = function(thre_tc) {
      if(!is.null(self$df)){
        if(min(self$df$tc_u_median) > thre_tc){
          return(FALSE)
        }else{
          return(TRUE)
        }
      }else{
        return(NA)
      }
    }
  )
)

#' @title Predictions class for curve values
#'
#' @description
#' Concrete class of Predictions for curve values.
#'
#' @name Predictions_Curve
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Predictions}}{}
#' }
#' 
#' @section Public fields:
#' 
#' @keywords predictions, Stan
#' @family Logistic4
Predictions_Curve <- R6::R6Class("Predictions_Curve",  inherit = Predictions,
                                 
  public = list(
    
    #' @description
    #' Set \code{y_pred_names}, \code{u_pred_names}, \code{df} and \code{xx} declared in the superclass
    #' using the argument \code{fit_info}.
    #' @param fit_info Fit_info class object.
    initialize = function(fit_info){
      if(!is.null(fit_info$fit)){
        self$y_pred_names <- paste0("y_pred[", 1:(fit_info$par_dims$y_pred),"]")
        self$u_pred_names <- paste0("u_pred[", 1:(fit_info$par_dims$u_pred),"]")
        super$set_predictions(fit_info$fit_summary, fit_info$input, fit_info$sf_tc)
      }
    }
  )
)

#' @title Predictions class for surface values
#'
#' @description
#' Concrete class of Predictions for surface values.
#'
#' @name Predictions_Surface
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Predictions}}{}
#' }
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{mx_surface_median}: }{matrix form of predicted u at 50 percentile.}
#'   \item{\code{mx_surface_lower}: }{matrix form of predicted u at 2.5 percentile.}
#'   \item{\code{mx_surface_upper}: }{matrix form of predicted u at 97.5 percentile.}
#'   \item{\code{mx_surface_median_tc}: }{matrix form of predicted T/C u at 50 percentile.}
#' }
#' @keywords predictions, Stan
#' @family BayeSyC4
Predictions_Surface <- R6::R6Class("Predictions_Surface",  inherit = Predictions,
                                   
  public = list(
    mx_surface_median = NULL,
    mx_surface_lower = NULL,
    mx_surface_upper = NULL,
    mx_surface_median_tc = NULL,
    
    #' @description
    #' Set \code{y_pred_names}, \code{u_pred_names} and \code{df} declared in the superclass
    #' using the argument \code{fit_info}.
    #' Set \code{mx_surface_median}, \code{mx_surface_lower}, \code{mx_surface_upper} and \code{mx_surface_median_tc}.
    #' @param fit_info Fit_info class object.
    initialize = function(fit_info){
      if(!is.null(fit_info$fit)){
        self$y_pred_names <- paste0("yc_pred[", 1:(fit_info$par_dims$yc_pred),"]")
        self$u_pred_names <- paste0("uc_pred[", 1:(fit_info$par_dims$uc_pred),"]")
        super$set_predictions(fit_info$fit_summary, fit_info$input, fit_info$sf_tc)
        # Matrix for surface values 
        self$mx_surface_median <- private$get_mx_surface(self$df$u_median, fit_info$input$v_idx1_new, fit_info$input$v_idx2_new)
        self$mx_surface_lower <- private$get_mx_surface(self$df$u_lower, fit_info$input$v_idx1_new, fit_info$input$v_idx2_new)
        self$mx_surface_upper <- private$get_mx_surface(self$df$u_upper, fit_info$input$v_idx1_new, fit_info$input$v_idx2_new)
        self$mx_surface_median_tc <- self$mx_surface_median * fit_info$sf_tc
      }
    }
  ),
  
  private = list(
    #' @description
    #' @param v_value vector of surface values 
    #' @param v_idx1 vector of first indices in sparse matrix.
    #' @param v_idx2 vector of second indices in sparse matrix
    #' @return matrix form of surface values.
    get_mx_surface = function(v_value, v_idx1, v_idx2){
      sp_surface <- Matrix::sparseMatrix(i=v_idx2, j=v_idx1, x=v_value) # i:y, j:x
      mx_surface <- as.matrix(sp_surface)
      return(mx_surface)
    }
  )
)
