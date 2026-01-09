
#' @title Data handler interface
#'
#' @description
#' Interface of data handling for Stan input. 
#'
#' @name Data_Handler
#' @docType class
#' 
#' @section Public fields:
#' 
#' @keywords data, Stan
Data_Handler <- R6::R6Class("Data_Handler",
  public = list(
    get_drc_input_data = function(assay_data, comp, params, metric) stop("I'm an abstract interface method."),
    get_bayesyc_input_data = function(assay_data, params, metric) stop("I'm an abstract interface method.")
  )
)



#' @title Data handler class 1
#'
#' @description
#' Concrete class of Data_Handler.
#'
#' @name Data_Handler_1
#' @docType class
#' 
#' @section Inherit:
#' \describe{
#'   \item{\code{Data_Handler}}{}
#' }
#' 
#' @section Public fields:
#' 
#' @keywords data, Stan
Data_Handler_1 <- R6::R6Class("Data_Handler_1",  inherit = Data_Handler,
                                 
  public = list(
    
    #' @description
    #' Generate and return an input list object for logistic curve.
    #' @param assay_data Assay_Data class object.
    #' @param comp target compound.
    #' @param params list object of parameters.
    #' @param metric target metric in multi-metric mode. NULL for single-metric mode.
    #' @return A list object.
    get_drc_input_data = function(assay_data, comp, params, metric){
      df <- assay_data$df_x[assay_data$df_x$compound == comp, ]
      df$conc <- as.numeric(df$conc)
      v_hc <- assay_data$df_hc[, metric]
      if(is.null(assay_data$df_lc)){
        v_lc <- NULL
      }else{
        v_lc <- assay_data$df_lc[, metric]
      }
        
      # Prediction points
      xmin <- min(df$conc)
      xmax <- max(df$conc)
      xx <- exp(seq(log(xmin), log(xmax), length.out=params$l_pred))
      # Input data
      ls_in_drc <- list(N=nrow(df), y=df$y, x=df$conc,
                        x_hat=median(df$conc), x_max=max(df$conc),
                        lb_log10_C=ifelse(is.null(params$lb_log10_C), -10, params$lb_log10_C),
                        ub_log10_C=ifelse(is.null(params$ub_log10_C), 6, params$ub_log10_C),
                        E0=median(assay_data$v_hc), hc_med=median(v_hc),
                        N_new=length(xx), x_new=xx, df=df,
                        N_hc=length(v_hc), hc=v_hc,
                        N_lc=length(v_lc), lc=v_lc)
      return(ls_in_drc)
    },
    
    #' Generate and return an input list object for BayeSyC.
    #' @param assay_data Assay_Data class object.
    #' @param params list object of parameters.
    #' @param metric target metric in multi-metric mode. NULL for single-metric mode.
    #' @return A list object.
    get_bayesyc_input_data = function(assay_data, params, metric){
      df1 <- assay_data$df_x1
      df2 <- assay_data$df_x2
      v_hc <- assay_data$df_hc[, metric]
      if(is.null(assay_data$df_lc)){
        v_lc <- NULL
      }else{
        v_lc <- assay_data$df_lc[, metric]
      }
      
      #------------------------------
      # Matrix-shape dose response curve
      if(params$use_zeros){
        dfc <- assay_data$df_all_matrix
        dfc[dfc$conc_1 == 0, "conc_1"] <- min(1e-8, 0.001 * min(dfc[dfc$conc_1 > 0, "conc_1"]))
        dfc[dfc$conc_2 == 0, "conc_2"] <- min(1e-8, 0.001 * min(dfc[dfc$conc_2 > 0, "conc_2"]))
      }else{
        #df_x_matrix <- assay_data$df_x_matrix
        #dfc <- df_x_matrix[(df_x_matrix$conc_1 !=0 & df_x_matrix$conc_2 !=0), ]
        dfc <- assay_data$df_x_matrix
      }
      
      # Prediction points
      df_x_matrix <- assay_data$df_x_matrix
      v_x1_new <- exp(seq(log(min(df_x_matrix$conc_1)), log(max(df_x_matrix$conc_1)), length.out=params$l_pred))
      v_x2_new <- exp(seq(log(min(df_x_matrix$conc_2)), log(max(df_x_matrix$conc_2)), length.out=params$l_pred))
      v_dc1_new <- c()
      v_dc2_new <- c()
      v_idx1_new <- c()
      v_idx2_new <- c()
      for(idx1 in seq(v_x1_new)){
        for(idx2 in seq(v_x2_new)){
          v_idx1_new <- c(v_idx1_new, idx1)
          v_idx2_new <- c(v_idx2_new, idx2)
          v_dc1_new <- c(v_dc1_new, v_x1_new[idx1])
          v_dc2_new <- c(v_dc2_new, v_x2_new[idx2])
        }
      }
      
      # Input data for BayeSyC and Dual dose response curve
      ls_in_bayesyc <- list(
        # For dual
        N1=nrow(df1), y1=df1[[metric]], x1=df1$conc,
        N2=nrow(df2), y2=df2[[metric]], x2=df2$conc,
        Nh=length(v_hc), hc=v_hc, hc_med=median(v_hc), hc_mean=mean(v_hc), 
        Nl=length(v_lc), lc=v_lc, 
        lb_log10_C=ifelse(is.null(params$lb_log10_C), -10, params$lb_log10_C), 
        ub_log10_C=ifelse(is.null(params$ub_log10_C), 6, params$ub_log10_C),
        E0=median(v_hc),
        # For combination
        Nc=nrow(dfc), yc=dfc[[metric]], dc1=dfc$conc_1, dc2=dfc$conc_2,
        r1r=ifelse(is.null(params$bayesyc_r1r), 0.0001, params$bayesyc_r1r), 
        r2r=ifelse(is.null(params$bayesyc_r2r), 0.0001, params$bayesyc_r2r), 
        gamma12=ifelse(is.null(params$bayesyc_gamma12), 1, params$bayesyc_gamma12),
        gamma21=ifelse(is.null(params$bayesyc_gamma21), 1, params$bayesyc_gamma21),
        Nc_new=length(v_dc1_new), dc1_new=v_dc1_new, dc2_new=v_dc2_new,
        v_idx1_new=v_idx1_new, v_idx2_new=v_idx2_new,  v_x1_new=v_x1_new, v_x2_new=v_x2_new,
        df_x_matrix=df_x_matrix,
        # Mainly for SynBa
        einf_beta_a=ifelse(is.null(params$bayesyc_einf_beta_a), 1, params$bayesyc_einf_beta_a), 
        einf_beta_b=ifelse(is.null(params$bayesyc_einf_beta_b), 1, params$bayesyc_einf_beta_b), 
        sigma_mu=ifelse(is.null(params$bayesyc_sigma_mu), 0, params$bayesyc_sigma_mu) 
      )
      return(ls_in_bayesyc)
      
    }
  )
)



#' @title Assay data class
#'
#' @description
#' Class of assay data providing functions for check, selection and making summary.
#'
#' @name Assay_Data
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{df_x}: }{data.frame of single-drug response data.}
#'   \item{\code{df_all_matrix}: }{data.frame of matrix-form observed data including points at a concentration of 0. NULL if the compound layout is not matrix-form.}
#'   \item{\code{df_x_matrix}: }{data.frame of matrix-form observed data. NULL if the compound layout is not matrix-form.}
#'   \item{\code{df_x1}: }{data.frame of single-drug response data.}
#'   \item{\code{df_x2}: }{data.frame of single-drug response data.}
#'   \item{\code{v_hc}: }{numeric vector of high control.}
#'   \item{\code{specimen}: }{specimen.}
#'   \item{\code{analysis}: }{analysis.}
#'   \item{\code{display_unit}: }{display_unit.}
#'   \item{\code{v_comp}: }{character vector of compounds.}
#'   \item{\code{v_metric}: }{character vector of metrics used.}
#'   \item{\code{data_handler}: }{Data_Handler class object that generates an input list for Stan sampling.}
#' }
#' @keywords assay, data
#' @family Data
#' 
#' @export
Assay_Data <- R6::R6Class("Assay_Data",

  public = list(
             df_x = NULL, 
             df_all_matrix = NULL,
             df_x_matrix = NULL,
             df_x1 = NULL, 
             df_x2 = NULL, 
             col_target = NULL,
             
             df_lc = NULL, df_hc = NULL,
             # Summary statistics
             v_hc_mean = NULL, v_hc_med = NULL, v_hc_sd = NULL, v_lc_mean = NULL, v_lc_med = NULL, v_lc_sd = NULL,
             
             specimen = NULL,
             analysis = NULL,
             v_comp = NULL,
             v_metric = NULL,
             data_handler = NULL,
             display_unit = NULL,
    
    #' @description
    #' Initialize fields by \code{data} that contains elements correspond to the fields.
    #' @param data list of meta data and observed data corresponding to one assay.
    initialize = function(data){
      self$specimen <- data$specimen
      self$analysis <- data$analysis
      self$df_x <- data$df_x
      self$df_all_matrix <- data$df_all_matrix
      self$df_x_matrix <- data$df_x_matrix
      self$df_x1 <- data$df_x1
      self$df_x2 <- data$df_x2
      self$df_lc <- data$df_lc
      self$df_hc <- data$df_hc
      self$v_comp = data$v_comp
      self$v_metric = data$v_metric
      self$data_handler <- data$data_handler
      self$display_unit <- data$analysis
    },
    
    #' @description
    #' Process data.
    #' @param select "target" for single target metric OR "multi" for multimetric.
    #' @param params list object of parameters.
    process_data = function(select, params){
      self$col_target <- params$col_target
      self$v_metric <- self$col_target
      return(NULL)
    },
    
    #' @description
    #' Generate and return an input list for Stan sampling of logistic curve.
    #' @param comp target compound.
    #' @param params list object of parameters.
    #' @param metric target metric in multi-metric mode. NULL for single-metric mode.
    #' @return A list object.
    get_drc_input_data = function(comp, params, metric){
      return(self$data_handler$get_drc_input_data(self, comp, params, metric))
    },
    
    #' @description
    #' Generate and return an input list for Stan sampling of BayeSyC.
    #' @param params list object of parameters.
    #' @param metric target metric in multi-metric mode. NULL for single-metric mode.
    #' @return A list object.
    get_bayesyc_input_data = function(params, metric){
      return(self$data_handler$get_bayesyc_input_data(self, params, metric))
    },
    
    #' @description
    #' Calculate summary statistics for high and low control values of \code{v_metric}.
    set_summary_data = function(){
      # Set mean, median and SD of HC and LC
      self$v_hc_mean <- apply(self$df_hc[, self$v_metric, drop = FALSE], 2, mean)
      self$v_hc_med <- apply(self$df_hc[, self$v_metric, drop = FALSE], 2, median)
      self$v_hc_sd <- apply(self$df_hc[, self$v_metric, drop = FALSE], 2, sd)
      self$v_lc_mean <- apply(self$df_lc[, self$v_metric, drop = FALSE], 2, mean)
      self$v_lc_med <- apply(self$df_lc[, self$v_metric, drop = FALSE], 2, median)
      self$v_lc_sd <- apply(self$df_lc[, self$v_metric, drop = FALSE], 2, sd)
    },
    
    #' @description
    #' Return meta data for stats output.
    #' @return A character vector of meta data.
    get_meta = function(){
      v_meta <- c(self$specimen)
      names(v_meta) <- c("Specimen")
      return(v_meta)
    },
    
    #' @description
    #' Pseudo function returning NULL for stats output.
    #' @return A character vector of judgement information.
    get_judgement = function(){
      return(NULL)
    },
    
    #' @description
    #' Pseudo function returning NULL for stats output.
    #' @param df data.frame of one particular drug response.
    #' @return A character vector of limit of detection information.
    get_lod_info = function(df){
      return(NULL)
    }
  )
)


