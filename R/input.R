
#' @title Class for loading single-drug dose response data
#'
#' @description
#' This class for manually loading single-drug dose response data
#'
#' @name Manual_Input_Single
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{ls_data}: }{list of Assay_Data class objects.}
#'   \item{\code{v_specimen}: }{character vector of specimen IDs.}
#'   \item{\code{v_analysis}: }{character vector of analysis IDs.}
#'   \item{\code{uniq_specimens}: }{character vector of unique specimen IDs.}
#'   \item{\code{uniq_compounds}: }{character vector of unique compounds.}
#'   \item{\code{n_uniq_specimens}: }{number of unique specimen IDs.}
#'   \item{\code{n_uniq_compounds}: }{number of unique compounds.}
#'   \item{\code{uniq_metrics}: }{character vector of unique metrics.}
#'   \item{\code{n_uniq_metrics}: }{number of unique metrics.}
#'   \item{\code{uniq_display_units}: }{character vector of unique display units.}
#' }
#' @keywords manual, single
#' @family Input
#'
#' @export
Manual_Input_Single <- R6::R6Class("Manual_Input_Single",
                                   
  public = list(
    ls_data = NULL,
    v_specimen = NULL,
    v_analysis = NULL,
    uniq_specimens = NULL,
    uniq_compounds = NULL,
    n_uniq_specimens = NULL,
    n_uniq_compounds = NULL,
    uniq_metrics = NULL,
    n_uniq_metrics = NULL,
    uniq_display_units = NULL,
    
    #' @description
    #' Load input data.frame and standardize data format for downstream analysis.
    #' @param params list object of parameters
    initialize = function(params, df_single, base_value=NULL){
      
      v_specimen <- unique(na.omit(df_single$specimen))
      v_compound <- NULL
      
      # Low control
      df_lc <- df_single %>% dplyr::filter(is.na(specimen))
      
      for(specimen in v_specimen){ 
        df <- df_single %>% dplyr::filter(specimen == specimen, !is.na(y))
        # High control
        df_hc <- df %>% dplyr::filter(conc == 0)
        if(nrow(df_hc) == 0){
          if(is.null(base_value)){
            df_hc <- df %>% dplyr::filter(conc == min(conc))
            warning(paste("No high control. Use values at minimum concentrations instead.\n"))
          }else{
            df_hc <- data.frame(y=base_value)
            warning(paste("No high control. Use specified base value instead:", base_value, "\n"))
          }
        }
        df_x <- df[df$conc != 0, ]
        v_comp <- unique(df_x$compound)
        
        # Store data
        assay_data <- Assay_Data$new(
          list(specimen=specimen, analysis=specimen,
               df_x=df_x, df_hc=df_hc, df_lc=df_lc, v_comp=v_comp,
               col_target=params$col_target,
               data_handler=Data_Handler_1$new())
        )
        
        self$ls_data <- c(self$ls_data, list(assay_data))
        self$v_specimen <- c(self$v_specimen, specimen)
        self$v_analysis <- c(self$v_analysis, specimen)
        v_compound <- c(v_compound, v_comp)
      }
      
      self$uniq_specimens <- unique(self$v_specimen)
      self$n_uniq_specimens <- length(self$v_specimen)
      self$uniq_compounds <- unique(v_compound)
      self$n_uniq_compounds <- length(self$uniq_compounds)
    },
    
    
    #' @description
    #' Process input data
    #' @param select what to be selected: "target" or "multi"
    #' @param params list object of parameters
    #' @return A character vector of error messages. NULL if no error.
    process_input_data = function(select, params){
      v_all_metric <- c()
      v_display_unit <- c()
      for(p in seq(self$ls_data)){
        self$ls_data[[p]]$process_data(select, params)
        v_all_metric <- c(v_all_metric, self$ls_data[[p]]$v_metric)
        v_display_unit <- c(v_display_unit, self$ls_data[[p]]$display_unit)
      }
      # All unique metrics across the whole input
      self$uniq_metrics <- unique(v_all_metric) 
      self$n_uniq_metrics <- length(self$uniq_metrics)
      # Display unit across the whole input
      self$uniq_display_units <- unique(v_display_unit)
      return(NULL)
    }
  )
)


#' @title Class for loading two-drug combination dose response data
#'
#' @description
#' This class for manually loading two-drug combination dose response data
#'
#' @name Manual_Input_Double
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{ls_data}: }{list of Assay_Data class objects.}
#'   \item{\code{v_specimen}: }{character vector of specimen IDs.}
#'   \item{\code{v_analysis}: }{character vector of analysis IDs.}
#'   \item{\code{uniq_display_units}: }{character vector of unique display units.}
#' }
#' @keywords manual, double
#' @family Input
#'
#' @export
Manual_Input_Double <- R6::R6Class("Manual_Input_Double",
                                   
  public = list(
    ls_data = NULL,
    v_specimen = NULL,
    v_analysis = NULL,
    
    #' @description
    #' Load input data.frame and standardize data format for downstream analysis.
    #' @param params list object of parameters
    initialize = function(params, df_matrix, base_value=NULL){
      
      v_specimen <- unique(df_matrix$specimen)
      for(sp in v_specimen){ 
        df_matrix_sp <- df_matrix %>% dplyr::filter(specimen == sp)
        # Low control
        df_lc <- df_matrix_sp %>% dplyr::filter(is.na(combination_name))
        # without Low control
        df_matrix_sp <- df_matrix_sp %>% dplyr::filter(!is.na(combination_name))
        v_set <- unique(df_matrix_sp$combination_name)
        for(set in v_set){ 
          analysis <- paste0(sp, " ", set)
          df <- df_matrix_sp %>% dplyr::filter(combination_name == set, !is.na(y))
          # Double
          df_x_matrix <- df %>% dplyr::filter(conc_1 > 0, conc_2 > 0)
          # Single
          compound_1 <- unique(df$compound_1); stopifnot(length(compound_1) == 1)
          compound_2 <- unique(df$compound_2); stopifnot(length(compound_2) == 1)
          df_x1 <- df[df$conc_1 != 0 & df$conc_2 == 0, ] %>% dplyr::rename(conc=conc_1, compound=compound_1)
          df_x2 <- df[df$conc_1 == 0 & df$conc_2 != 0, ] %>% dplyr::rename(conc=conc_2, compound=compound_2)
          # High control
          df_hc <- df %>% dplyr::filter(conc_1 == 0, conc_2 == 0)
          if(nrow(df_hc) == 0){
            if(is.null(base_value)){
              df_hc <- df %>% dplyr::filter(conc_1 == min(conc_1), conc_2 == min(conc_2))
              warning(paste("No high control. Use values at minimum concentrations instead.\n"))
            }else{
              df_hc <- data.frame(y=base_value)
              warning(paste("No high control. Use specified base value instead:", base_value, "\n"))
            }
          }
          
          # Store data
          assay_data <- Assay_Data$new(
            list(specimen=sp, analysis=analysis,
                 df_x1=df_x1, df_x2=df_x2, df_x_matrix=df_x_matrix, df_all_matrix=df,
                 df_hc=df_hc, df_lc=df_lc,
                 data_handler=Data_Handler_1$new())
          )
            
          self$ls_data <- c(self$ls_data, list(assay_data))
          self$v_specimen <- c(self$v_specimen, sp)
          self$v_analysis <- c(self$v_analysis, analysis)
        }
      }
    },
    
    
    #' @description
    #' Process input data
    #' @param select what to be selected: "target" or "multi"
    #' @param params list object of parameters
    #' @return A character vector of error messages. NULL if no error.
    process_input_data = function(select, params){
      for(p in seq(self$ls_data)){
        self$ls_data[[p]]$process_data(select, params)
      }
      return(NULL)
    }
  )
)



#' @title Class for loading ONeil data
#'
#' @description
#' This class handle O'Neil et al. 2016 data.
#'
#' @name Whole_Input_ONeil
#' @docType class
#'
#' @section Public fields:
#' \describe{
#'   \item{\code{ls_data}: }{list of Assay_Data class objects.}
#'   \item{\code{v_specimen}: }{character vector of specimen IDs.}
#'   \item{\code{v_analysis}: }{character vector of analysis IDs.}
#'   \item{\code{n_uniq_metrics}: }{number of unique metrics.}
#' }
#' @keywords ONeil, double
#' @family Input
#'
#' @export
Whole_Input_ONeil <- R6::R6Class("Whole_Input_ONeil",
                                   
  public = list(
    ls_data = NULL,
    v_specimen = NULL,
    v_analysis = NULL,
    
    #' @description
    #' Load input data.frame and standardize data format for downstream analysis.
    #' @param params list object of parameters
    initialize = function(params, df_double, df_single){
      # Double
      df_double[, private$cols_viability_4] <- sapply(df_double[, private$cols_viability_4], as.numeric) # NULL will be NA here
      if(grepl("gamma", params$model_bayesyc)){
        # For gamma likelihood function
        # Set viability = 0 if viability < 0
        df_double <- df_double %>% dplyr::mutate(across(all_of(private$cols_viability_4), ~ pmax(., 1e-4))) #ifelse(. <= 0, 1e-5, .)))
      }
      # Long format
      dfl_double <- tidyr::pivot_longer(df_double, cols=private$cols_viability_4, names_to="replicate", values_to="y" )
      dfl_double <- dfl_double[!is.na(dfl_double$y), ] # excluding NAs
      # rename
      lookup <- c(compound_1="drugA_name", compound_2="drugB_name", conc_1="drugA_conc_uM", conc_2="drugB_conc_uM")
      dfl_double <- dplyr::rename(dfl_double, dplyr::all_of(lookup))
      
      # Single
      df_single[, private$cols_viability_6] <- sapply(df_single[, private$cols_viability_6], as.numeric) # NULL will be NA here
      if(grepl("gamma", params$model_bayesyc)){
        # For gamma likelihood function
        # Set viability = 0 if viability < 0
        df_single <- df_single %>% dplyr::mutate(across(all_of(private$cols_viability_6), ~ pmax(., 1e-4))) #ifelse(. <= 0, 1e-5, .)))
      }
      # Long format
      dfl_single <- tidyr::pivot_longer(df_single, cols=private$cols_viability_6, names_to="replicate", values_to="y" )
      dfl_single <- dfl_single[!is.na(dfl_single$y), ] # excluding NAs
      colnames(dfl_single)[c("drug_name", "drug_conc_uM")] <- c("compound", "conc")
      # rename
      lookup <- c(compound="drug_name", conc="drug_conc_uM")
      dfl_single <- dplyr::rename(dfl_single, dplyr::all_of(lookup))
      
      v_cell_line <- unique(df_double$cell_line)
      for(cell_line in v_cell_line){ 
        v_set <- unique(dfl_double[dfl_double$cell_line == cell_line, "combination_name"])[[1]]
        for(set in v_set){ 
          analysis <- paste0(cell_line, " ", set)
          df_x_matrix <- dfl_double[dfl_double$cell_line == cell_line & dfl_double$combination_name == set, ]
          compound_1 <- unique(df_x_matrix$compound_1); stopifnot(length(compound_1) == 1)
          compound_2 <- unique(df_x_matrix$compound_2); stopifnot(length(compound_2) == 1)
          df_x1 <- dfl_single[dfl_single$cell_line == cell_line & dfl_single$compound == compound_1, ]
          df_x2 <- dfl_single[dfl_single$cell_line == cell_line & dfl_single$compound == compound_2, ]
          n_x1 <- nrow(df_x1); n_x2 <- nrow(df_x2)
          if(nrow(df_x1) == 0 | nrow(df_x2) == 0) next
          # df_all_matrix
          df_all_matrix <- rbind(df_x_matrix,
                                 data.frame(BatchID=df_x1$BatchID, cell_line=rep(cell_line, n_x1),
                                            compound_1=rep(compound_1, n_x1), conc_1=df_x1$conc,
                                            compound_2=rep(compound_2, n_x1), conc_2=rep(0, n_x1),
                                            combination_name=rep(set, n_x1), `mu/muMax`=df_x1$`mu/muMax`,
                                            `X/X0`=df_x1$`X/X0`, replicate=df_x1$replicate, y=df_x1$y, check.names=F),
                                 data.frame(BatchID=df_x2$BatchID, cell_line=rep(cell_line, n_x2),
                                            compound_1=rep(compound_1, n_x2), conc_1=rep(0, n_x2),
                                            compound_2=rep(compound_2, n_x2), conc_2=df_x2$conc,
                                            combination_name=rep(set, n_x2), `mu/muMax`=df_x2$`mu/muMax`,
                                            `X/X0`=df_x2$`X/X0`, replicate=df_x2$replicate, y=df_x2$y, check.names=F),
                                 data.frame(BatchID=0, cell_line=cell_line,
                                            compound_1=compound_1, conc_1=0, compound_2=compound_2, conc_2=0,
                                            combination_name=set, `mu/muMax`=1, `X/X0`=1, replicate="viability0", y=1, check.names=F)
                                )
          
          # Store data
          assay_data <- Assay_Data$new(
            list(specimen=cell_line, analysis=analysis,
                 df_x1=df_x1, df_x2=df_x2, df_x_matrix=df_x_matrix, df_all_matrix=df_all_matrix, 
                 v_hc=1, df_hc=data.frame(y=1), 
                 v_comp=c(compound_1, compound_2),
                 data_handler=Data_Handler_1$new())
          )
          
          self$ls_data <- c(self$ls_data, list(assay_data))
          self$v_specimen <- c(self$v_specimen, cell_line)
          self$v_analysis <- c(self$v_specimen, analysis)
        }
      }
    },
    
    
    #' @description
    #' Process input data
    #' @param select what to be selected: "target" or "multi"
    #' @param params list object of parameters
    #' @return A character vector of error messages. NULL if no error.
    process_input_data = function(select, params){
      for(p in seq(self$ls_data)){
        self$ls_data[[p]]$process_data(select, params)
      }
      return(NULL)
    }
  ),
  
  private = list(
    cols_viability_4 = c("viability1", "viability2", "viability3", "viability4"),
    cols_viability_6 = c("viability1", "viability2", "viability3", "viability4", "viability5", "viability6")
  )
)



