
#' @title Class for all fitting results
#'
#' @description
#' Class for all fitting results. This class should be inherited by concrete class.
#'
#' @name Whole_Output
#' @docType class
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{ls_data}: }{list of Plate_Data class objects.}
#'   \item{\code{params}: }{knit param list.}
#'   \item{\code{v_n_fit}: }{character of plate ID.}
#' }
#' @keywords fit, Stan
#' @family Output
#' 
#' @export
Whole_Output <- R6::R6Class("Whole_Output", 

  public = list(
    ls_data = NULL,
    params = NULL,
    v_n_fit = NULL,
    
    #' @description
    #' Execute MCMC sampling. If the result does not converge, it re-runs \code{self$params$n_retry} times.
    #' @param stan_model model (stanmodel or CmdStanModel) object.
    #' @param ls_in list object of input data.
    #' @param id ID of sampling
    #' @param mcmc_timeout max CPU time (time used for all cores are summed up) for MCMC sampling
    execute_stan_sampling = function(stan_model, ls_in, id, mcmc_timeout=20000){
      set.seed(777)
      for(n_run in seq(1 + self$params$n_retry)){
        #warning(paste(title, n_run))
        warning(paste(id, n_run))
        v_Rhat <- 99
        if(class(stan_model)[1] == "CmdStanModel"){
          if(is.null(self$params$cmdstanr_tmp)){
            tmp_output_dir <- NULL
          }else{
            # CmdStanR temporary directory
            # to avoid error from internal read_cmdstan_csv()
            tmp_output_dir <- file.path(self$params$cmdstanr_tmp, id)
            if(!dir.exists(tmp_output_dir)) dir.create(tmp_output_dir, recursive=TRUE)
          }
          # Perform CmdStanR sampling
          start_time <- proc.time()
          fit <- try(R.utils::withTimeout(
                      stan_model$sample(data=ls_in[!sapply(ls_in, is.null)], # remove NULL items
                      chains=self$params$n_chains, 
                      parallel_chains=self$params$n_cores, 
                      iter_warmup=self$params$n_iter, 
                      iter_sampling=self$params$n_iter,
                      refresh=self$params$refresh,
                      output_dir=tmp_output_dir),
                      timeout=mcmc_timeout, cpu=mcmc_timeout
                    ))
          end_time <- proc.time()
          elapsed_time <- end_time - start_time
          
          if(!is.null(fit)){
            if(inherits(fit, "try-error")){
              cat(paste0(id, ": Error or time out."), "\n")
              cat(paste0("Elapsed time: ", signif(elapsed_time["elapsed"])), "\n")
              fit <- NULL
            }else if(all(fit$return_codes() == 0)){
              v_Rhat <- fit$summary()[,"rhat"]
            }else{
              cat(paste0(id, ": Failed."), "\n")
              rm(fit)
              fit <- NULL
            }
          }
          # Remove temporary CmdStanR directory
          if(!is.null(tmp_output_dir)){
            if(dir.exists(tmp_output_dir)) unlink(tmp_output_dir, recursive = TRUE)
          }
        }else{
          # Perform RStan sampling
          fit <- rstan::sampling(stan_model, data=ls_in[!sapply(ls_in, is.null)],
                                 chains=self$params$n_chains, 
                                 cores=self$params$n_cores, 
                                 iter=self$params$n_iter,
                                 refresh=self$params$refresh)
          if(length(fit@sim) > 0){
            v_Rhat <- rstan::summary(fit)$summary[,"Rhat"]
          }else{
            fit <- NULL
          }
        }
        cat(paste("  max Rhat:", max(v_Rhat, na.rm=T)), "\n")
        # Finish sampling if all Rhats are OK.
        if(all(v_Rhat < self$params$thre_rhat, na.rm=T)) break
        # Finish sampling if sampling fails.
        if(is.null(fit)) break
      }
      return(fit)
    }
  )
)



#' @title Class for all multi-metric logistic curve fitting results 
#'
#' @description
#' This class conduct multi-metric logistic curve fitting and processing of the results for output.
#'
#' @name Whole_Output_Multi
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Whole_Output}}{}
#' }
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{stan_model_drc}: }{model (stanmodel or CmdStanModel) object of mono-phasic logistic curve.}
#'   \item{\code{stan_model_biphasic}: }{model (stanmodel or CmdStanModel) object of bi-phasic logistic curve.}
#'   \item{\code{ls_fit_multi}: }{list of Fit_Info class objects.}
#'   \item{\code{ls_res_multi}: }{list of processed fitting results.}
#'   \item{\code{df_stats_multi}: }{data.frame of stats.}
#'   \item{\code{ls_out_multi}: }{list of output-ready data for plots.}
#' }
#' 
#' @keywords fit, Stan
#' @family Output
#' 
#' @export
Whole_Output_Multi <- R6::R6Class("Whole_Output_Multi", inherit = Whole_Output,
                                  
  public = list(
    stan_model_drc = NULL,
    stan_model_biphasic = NULL,
    ls_fit_multi = NULL,
    ls_res_multi = NULL,
    df_stats_multi = NULL,
    ls_out_multi = NULL,
    
    #' @description
    #' Run MCMC sampling and process the fitting results for output.
    #' If argument \code{$run_all} FALSE, just making an instance.
    #' @param whole_input Whole_Input class object.
    #' @param stan_model_drc model (stanmodel or CmdStanModel) object of mono-phasic logistic curve.
    #' @param stan_model_biphasic  model (stanmodel or CmdStanModel) object of bi-phasic logistic curve.
    #' @param params list object of parameters
    #' @param run_all TRUE (default). FALSE for not running analyses and processes.
    #' @param bi_first FLASE (default). TRUE to prioritize bi-phasic logistic curve.
    initialize = function(whole_input, stan_model_drc, stan_model_biphasic, params, run_all=TRUE, bi_first=FALSE){
      
      self$ls_data <- whole_input$ls_data
      self$params <- params
      self$stan_model_drc <- stan_model_drc
      self$stan_model_biphasic <- stan_model_biphasic

      if(run_all){
        self$perform_fitting_multi(bi_first)
      
        if(all(self$v_n_fit > 0)){
          self$process_fit_multi()
          self$make_stats_multi()
          self$make_plot_data_multi()
        }
      }
    },

    #' @description
    #' Run MCMC sampling and hold the fitting results as \code{self$ls_fit_multi}.    
    #' @param bi_first FLASE (default). TRUE to prioritize bi-phasic logistic curve.
    perform_fitting_multi = function(bi_first=FALSE){
      
      for(p in seq(self$ls_data)){
        
        ls_f_comp <- NULL
        n_fit <- 0
        # Loop for each compound
        for(c in seq(self$ls_data[[p]]$v_comp)){
          comp <- self$ls_data[[p]]$v_comp[c]
          
          ls_f_metric <- NULL
          for(m in seq(self$ls_data[[p]]$v_metric)){
            metric <- self$ls_data[[p]]$v_metric[m]
            # Input data
            ls_in <- self$ls_data[[p]]$get_drc_input_data(comp, self$params, metric)
            
            if(ls_in$N > 1){
              #--------------------
              # Monophasic
              title_1 <- paste(self$ls_data[[p]]$analysis, comp, metric, "Monophasic Sampling", sep=", ")
              # Perform sampling
              id <- paste("drc", p, c, m, sep="_")
              cat(paste("Start:", id, title_1), "\n")
              fit1 <- super$execute_stan_sampling(self$stan_model_drc, ls_in, id)
              cat(paste("Finish:", id, title_1), "\n")
              # Create fit info object
              fit_drc1 <- Fit_MonoPhasic4$new(fit=fit1, input=ls_in, model=self$stan_model_drc)
              
              if(!is.null(self$stan_model_biphasic)){
                #--------------------
                # Biphasic
                title_2 <- paste(self$ls_data[[p]]$analysis, comp, metric, "Biphasic Sampling", sep=", ")
                # Perform sampling
                id <- paste("biphasic", p, c, m, sep="_")
                cat(paste("Start:", id, title_2), "\n")
                fit2 <- super$execute_stan_sampling(self$stan_model_biphasic, ls_in, id)
                cat(paste("Finish:", id, title_2), "\n")
                # Create fit info object
                fit_drc2 <- Fit_BiPhasic4$new(fit=fit2, input=ls_in, model=self$stan_model_biphasic)
                
                # Select model
                fit_drc <- fit_drc1$select_fit_info(fit_drc2, self$params$thre_rhat, bi_first) # check_fit is run internally
              }else{
                fit_drc <- fit_drc1
              }
              
            }else{
              warning(paste(self$ls_data[[p]]$analysis, comp, metric, ": Not enough input data.\n"))
              fit_drc <- Fit_MonoPhasic4$new(fit=NULL, input=ls_in, model=self$stan_model_drc)
            }
            
            # Check fit
            fit_drc$check_fit(self$params$thre_rhat)
            
            # Store in list
            ls_f_metric <- c(ls_f_metric, list(list(fit_drc=fit_drc)))
            # Number of fits
            if(!is.null(fit_drc$fit)) n_fit <- n_fit + 1
          }
          names(ls_f_metric) <- self$ls_data[[p]]$v_metric
          ls_f_comp <- c(ls_f_comp, list(ls_f_metric))
        }
        names(ls_f_comp) <- self$ls_data[[p]]$v_comp
        self$ls_fit_multi <- c(self$ls_fit_multi, list(ls_f_comp))
        self$v_n_fit <- c(self$v_n_fit, n_fit)
      }
    },

    #' @description
    #' Process the fitting results and hold processed results as \code{self$ls_res_multi}.
    process_fit_multi = function(){
      
      for(p in seq(self$ls_data)){
        
        ls_r_comp <- NULL
        for(comp in self$ls_data[[p]]$v_comp){
          
          ls_r_metric <- NULL
          for(metric in self$ls_data[[p]]$v_metric){
            # Stan fit info object
            fit_drc <- self$ls_fit_multi[[p]][[comp]][[metric]]$fit_drc
            # Data
            df <- fit_drc$input$df
            
            # Set predictions
            fit_drc$set_predictions()
            # Check minimum T/C
            fit_drc$check_tc(self$params$thre_tc)
            # Check IC50
            fit_drc$check_ic50(self$params$thre_ic50)
            
            # Set estimated pi (probability of outlier)
            df$tc <- df[, metric] * fit_drc$sf_tc
            if(!is.null(fit_drc$get_pi_names())){
              df$pi <- fit_drc$estimates$df[fit_drc$get_pi_names(), "50%"]
              df$well_info <- paste(df$well_info, signif(df$pi, 2))
            }else{
              df$pi <- rep(NA, nrow(df))
            }
            
            # Store in list
            ls_this <- list(p=p, file_no=paste0("No",as.character(p)), specimen=self$ls_data[[p]]$specimen, 
                            comp=comp, metric=metric, data=df, fit_drc=fit_drc)
            ls_r_metric <- c(ls_r_metric, list(ls_this))
          }
          names(ls_r_metric) <- self$ls_data[[p]]$v_metric
          ls_r_comp <- c(ls_r_comp, list(ls_r_metric))
        }
        names(ls_r_comp) <- self$ls_data[[p]]$v_comp
        self$ls_res_multi <- c(self$ls_res_multi, list(ls_r_comp))
      }
    },
    
    #' @description
    #' Make stats from the processed results and hold them as \code{self$df_stats_multi}.
    make_stats_multi = function(){
      
      ls_stats <- NULL
      v_name <- NULL
      for(p in seq(self$ls_res_multi)){
        for(c in seq(self$ls_res_multi[[p]])){
          for(metric in self$ls_data[[p]]$v_metric){
            
            fit_drc <- self$ls_res_multi[[p]][[c]][[metric]]$fit_drc
            comp <- self$ls_res_multi[[p]][[c]][[metric]]$comp
            df <- fit_drc$input$df
            
            # Stats for each plate - compound - metric
            v_comp_metric <- c(comp, metric)
            names(v_comp_metric) <- c("Compound", "Metric")
            
            # Stats to output
            v_stats <- c(self$ls_data[[p]]$get_meta(), v_comp_metric, self$ls_data[[p]]$get_judgement(),
                         self$ls_data[[p]]$get_lod_info(df), fit_drc$get_judgement(), fit_drc$stats$values)
      
            # Common names 
            if(is.null(v_name)){
              v_name <- names(v_stats)
            }else{
              v_name <- intersect(v_name, names(v_stats))
            }
            ls_stats <- c(ls_stats, list(v_stats))
          }
        }
      }
      
      # Update to use common names
      for(l in seq(ls_stats)){
        ls_stats[[l]] <- ls_stats[[l]][v_name]
      }
      self$df_stats_multi <- as.data.frame(t(as.data.frame(ls_stats, col.names=seq(length(ls_stats)), 
                                                           check.names=F)), check.names=F)
    },
    
    #' @description
    #' Make plot data from the processed results and hold them as \code{self$ls_out_multi}.
    make_plot_data_multi = function(){
      
      for(p in seq(self$ls_res_multi)){
        for(c in seq(self$ls_res_multi[[p]])){
          df_obs <- NULL
          df_pred <- NULL
          df_ic50 <- NULL
          v_metric <- names(self$ls_res_multi[[p]][[c]])
          
          for(metric in v_metric){
            fit_drc <- self$ls_res_multi[[p]][[c]][[metric]]$fit_drc
            estimates <- fit_drc$estimates
            predictions <- fit_drc$predictions
            
            data <- self$ls_res_multi[[p]][[c]][[metric]]$data
            comp <- self$ls_res_multi[[p]][[c]][[metric]]$comp
            
            #------------------------------
            # Observed values
            df_o <- data.frame(metric=rep(metric, nrow(data)), well=data$well, well_info=data$well_info, 
                               conc=data$conc, y=data[, metric], tc=data$tc)
            
            if(is.null(df_obs)){ df_obs <- df_o }else{ df_obs <- dplyr::bind_rows(df_obs, df_o) }
            
            #------------------------------
            # Predicted values
            if(fit_drc$success_fit){
              df_p <- transform(predictions$df, xx=predictions$xx,
                                metric=rep(metric, length(predictions$xx)), check.names=F)
              
              if(is.null(df_pred)){ df_pred <- df_p }else{ df_pred <- dplyr::bind_rows(df_pred, df_p) }
            }
            
            #------------------------------
            # Estimated IC50
            if(fit_drc$success_fit & fit_drc$ok_IC50){
              df_i <- transform(estimates$df["C", ], 
                                metric=metric, text="", check.names=F)
            }else{
              df_i <- data.frame(`2.5%`=NA, `50%`=0, `97.5%`=NA, 
                                 metric=metric, text="N.D.", check.names=F)
            }
            if(is.null(df_ic50)){ df_ic50 <- df_i }else{ df_ic50 <- dplyr::bind_rows(df_ic50, df_i) }
            
          }
          # n_plate x n_compound different data.frames are stored 
          # n_plate x n_compound figures will be displayed
          ls_this <- list(df_obs=df_obs, df_pred=df_pred, df_ic50=df_ic50,
                          no=as.character(p), analysis=self$ls_data[[p]]$analysis,
                          specimen=self$ls_data[[p]]$specimen, compound=comp)
          self$ls_out_multi <- c(self$ls_out_multi, list(ls_this))
        }
      }
    }
  )
)
  

#' @title Class for all BayeSyC fitting results 
#'
#' @description
#' This class conduct BayeSyC fitting and processing of the results for output.
#'
#' @name Whole_Output_BayeSyC
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Whole_Output}}{}
#' }
#' 
#' @section Public fields:
#' \describe{
#'   \item{\code{stan_model_drc}: }{model (stanmodel or CmdStanModel) object of mono-phasic logistic curve.}
#'   \item{\code{stan_model_biphasic}: }{model (stanmodel or CmdStanModel) object of bi-phasic logistic curve.}
#'   \item{\code{ls_fit_bayesyc}: }{list of Fit_Info class objects.}
#'   \item{\code{ls_res_bayesyc}: }{list of processed fitting results.}
#'   \item{\code{df_stats_bayesyc}: }{data.frame of stats.}
#' }
#' 
#' @keywords fit, Stan
#' @family Output
#' 
#' @export
Whole_Output_BayeSyC <- R6::R6Class("Whole_Output_BayeSyC", inherit = Whole_Output,
                                  
  public = list(
    stan_model_dual = NULL,
    stan_model_bayesyc = NULL,
    ls_fit_bayesyc = NULL,
    ls_res_bayesyc = NULL,
    df_stats_bayesyc = NULL,
    
    #' @description
    #' Run sampling and process the fitting results for output.
    #' If argument \code{$run_all} FALSE, just making an instance.
    #' @param whole_input Whole_Input class object.
    #' @param stan_model_dual model (stanmodel or CmdStanModel) object of dual mono-phasic logistic curves.
    #' @param stan_model_bayesyc  model (stanmodel or CmdStanModel) object of BayeSyC response surface.
    #' @param params list object of parameters
    #' @param run_all TRUE (default). FALSE for not running analyses and processes.
    initialize = function(whole_input, stan_model_dual, stan_model_bayesyc, params, run_all=TRUE){
      
      self$ls_data <- whole_input$ls_data
      self$params <- params
      self$stan_model_dual <- stan_model_dual
      self$stan_model_bayesyc <- stan_model_bayesyc
      
      if(run_all){
        self$perform_fitting_bayesyc()
        
        if(all(self$v_n_fit > 0)){
          self$process_fit_bayesyc()
          self$make_stats_bayesyc()
        }
      }
    },

    #' @description
    #' Run perform_fitting_multi_bayesyc_1 or perform_fitting_multi_bayesyc_2
    perform_fitting_bayesyc = function(){
      if(is.null(self$stan_model_dual)){
        self$perform_fitting_bayesyc_1()
      }else{
        self$perform_fitting_bayesyc_2()
      }
    },
    
    #' @description
    #' Run MCMC sampling and hold the fitting results as \code{self$ls_fit_multi_bayesyc}.
    #' Case of no Dual.
    perform_fitting_bayesyc_1 = function(){
      
      for(p in seq(self$ls_data)){
        
        ls_f_metric <- NULL
        for(m in seq(self$ls_data[[p]]$v_metric)){
          metric <- self$ls_data[[p]]$v_metric[m]
          
          # Input data
          ls_in <- self$ls_data[[p]]$get_bayesyc_input_data(self$params, metric)
          
          title_m <- paste(self$ls_data[[p]]$analysis, metric, "BayeSyC Sampling", sep=", ")
          # Perform sampling
          id <- paste("bayesyc", p, m, sep="_")
          cat(paste("Start:", id, title_m), "\n")
          fit_m <- super$execute_stan_sampling(self$stan_model_bayesyc, ls_in, id)
          cat(paste("Finish:", id, title_m), "\n")
          # Create fit info object
          fit_bayesyc <- Fit_BayeSyC4$new(fit=fit_m, input=ls_in, model=self$stan_model_bayesyc)
          # Check fit
          fit_bayesyc$check_fit(self$params$thre_rhat)
          
          # Store in list
          ls_f_metric <- c(ls_f_metric, list(list(fit_bayesyc=fit_bayesyc)))
        }
        names(ls_f_metric) <- self$ls_data[[p]]$v_metric
        self$ls_fit_bayesyc <- c(self$ls_fit_bayesyc, list(ls_f_metric))
      }
    },
    
    
    #' @description
    #' Run MCMC sampling and hold the fitting results as \code{self$ls_fit_multi_bayesyc}.
    #' Case of with Dual followed by BayeSyC.
    perform_fitting_bayesyc_2 = function(){
      
      for(p in seq(self$ls_data)){
        
        ls_f_metric <- NULL
        for(m in seq(self$ls_data[[p]]$v_metric)){
          metric <- self$ls_data[[p]]$v_metric[m]
          
          #------------------------------
          # Dual
          ls_in <- self$ls_data[[p]]$get_bayesyc_input_data(self$params, metric)
          title_d <- paste(self$ls_data[[p]]$analysis, metric, "Dual Sampling", sep=", ")
          # Perform sampling for Dual
          id <- paste("dual", p, m, sep="_")
          cat(paste("Start:", id, title_d), "\n")
          fit_d <- super$execute_stan_sampling(self$stan_model_dual, ls_in, id)
          cat(paste("Finish:", id, title_d), "\n")
          # Create fit info object
          fit_dual <- Fit_Dual$new(fit=fit_d, input=ls_in, model=self$stan_model_dual)
          # Check fit
          fit_dual$check_fit(self$params$thre_rhat)
        
          if(!fit_dual$success_fit){
            warning(paste0(title_d, ": Fitting failed (does not satisfy the criteria).\n"))
            # Create NULL fit info object
            fit_bayesyc <- Fit_BayeSyC4$new(fit=NULL, input=ls_in, model=self$stan_model_bayesyc)
            # Check fit
            fit_bayesyc$check_fit(self$params$thre_rhat)
          }else{
            #------------------------------
            # BayeSyC
            # Input data for BayeSyC using parameters from Posterior distributions of Dual
            ls_in_bayesyc <- c(ls_in, fit_dual$get_bayesyc_params())
          
            title_m <- paste(self$ls_data[[p]]$analysis, metric, "BayeSyC Sampling", sep=", ")
            # Perform sampling for BayeSyC
            id <- paste("bayesyc", p, m, sep="_")
            cat(paste("Start:", id, title_m), "\n")
            fit_m <- super$execute_stan_sampling(self$stan_model_bayesyc, ls_in_bayesyc, id)
            cat(paste("Finish:", id, title_m), "\n")
            # Create fit info object
            fit_bayesyc <- Fit_BayeSyC4$new(fit=fit_m, input=ls_in_bayesyc, model=self$stan_model_bayesyc)
            # Check fit
            fit_bayesyc$check_fit(self$params$thre_rhat)
          }
        
          # Store in list
          ls_f_metric <- c(ls_f_metric, list(list(fit_bayesyc=fit_bayesyc, fit_dual=fit_dual)))
        }
        names(ls_f_metric) <- self$ls_data[[p]]$v_metric
        self$ls_fit_bayesyc <- c(self$ls_fit_bayesyc, list(ls_f_metric))
      }
    },
    
    #' @description
    #' Process the fitting results and hold processed results as \code{self$ls_res_multi}.
    process_fit_bayesyc = function(){
      
      for(p in seq(self$ls_data)){
        
        ls_r_metric <- NULL
        for(metric in self$ls_data[[p]]$v_metric){
          
          # Stan fit info object
          fit_bayesyc <- self$ls_fit_bayesyc[[p]][[metric]]$fit_bayesyc
          # Data
          df_x_matrix <- fit_bayesyc$input$df_x_matrix
          v_comb_comp <- c(unique(df_x_matrix$compound_1), unique(df_x_matrix$compound_2))
          
          # Set predictions
          fit_bayesyc$set_predictions()
          # Check minimum T/C
          fit_bayesyc$check_tc(self$params$thre_tc)
          # Set T/C of observed data
          df_x_matrix$tc <- df_x_matrix$y * fit_bayesyc$sf_tc
          
          # Store in list
          ls_this <- list(p=p, file_no=paste0("No",as.character(p)), specimen=self$ls_data[[p]]$specimen, 
                          df_x_matrix=df_x_matrix, comp=v_comb_comp, metric=metric, fit_bayesyc=fit_bayesyc)
          
          ls_r_metric <- c(ls_r_metric, list(ls_this))
        }
        names(ls_r_metric) <- self$ls_data[[p]]$v_metric
        self$ls_res_bayesyc <- c(self$ls_res_bayesyc, list(ls_r_metric))
      }
    },
    
    #' @description
    #' Make stats from the processed results and hold them as \code{self$df_stats_multi}.
    make_stats_bayesyc = function(){
      
      ls_stats <- NULL
      v_name <- NULL
      for(p in seq(self$ls_data)){
        for(metric in self$ls_data[[p]]$v_metric){
          
          fit_bayesyc <- self$ls_res_bayesyc[[p]][[metric]]$fit_bayesyc
          v_comb_comp <- self$ls_res_bayesyc[[p]][[metric]]$comp
          df <- self$ls_data[[p]]$df_x
        
          # Stats for each plate
          v_comp_metric <- c(v_comb_comp[1], v_comb_comp[2], metric)
          names(v_comp_metric) <- c("Compound_1", "Compound_2", "Metric")
        
          # Stats to output
          v_stats <- c(self$ls_data[[p]]$get_meta(), v_comp_metric, self$ls_data[[p]]$get_judgement(), 
                       self$ls_data[[p]]$get_lod_info(df), fit_bayesyc$get_judgement(), fit_bayesyc$stats$values)
        
          # Common names 
          if(is.null(v_name)){
            v_name <- names(v_stats)
          }else{
            v_name <- intersect(v_name, names(v_stats))
          }
          ls_stats <- c(ls_stats, list(v_stats))
        }
      }
        # Update to use common names
      for(l in seq(ls_stats)){
        ls_stats[[l]] <- ls_stats[[l]][v_name]
      }
      self$df_stats_bayesyc <- as.data.frame(t(as.data.frame(ls_stats, col.names=seq(length(ls_stats)), 
                                                           check.names=F)), check.names=F)
    }
  )
)
  
