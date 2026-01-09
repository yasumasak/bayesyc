
#' Get legend value for second key
#'
#' @param df data.frame used for plot
#' @param v_def_value definition of value used for manual setting
#' @param first_key first key
#' @param second_key second key
#'
#' @return v_value legend value for second key
#'
#' @export
get_second_legend_value <- function(df, v_def_value, first_key, second_key){
  mx_key <- unique(df[, c(first_key, second_key)])
  v_value <- c()
  pre_key <- ""
  idx <- 1

  for(i in seq(nrow(mx_key))){
    key = mx_key[i, first_key]
    if(pre_key != key){
      # If first key changes, the first value is used.
      idx <- 1
    }else{
      # If first key does not change (second key changes),
      # the next value is used.
      idx <- idx + 1
    }
    v_value <- c(v_value, v_def_value[idx])
    pre_key <- key
  }
  return(v_value)
}


#' Modify ggplotly legend
#'
#' Remove unnecessary parentheses "(*****,1)" and put space after comma
#'
#' @param obj_ggplotly ggplotly object
#'
#' @return obj_ggplotly ggplotly object
#'
#' @export
modify_ggplotly_legend <- function(obj_ggplotly){
  for (i in seq(obj_ggplotly$x$data)){
    if (!is.null(obj_ggplotly$x$data[[i]]$name)){
      obj_ggplotly$x$data[[i]]$name = gsub('^\\(|,+\\d+\\)$', '', obj_ggplotly$x$data[[i]]$name)
      obj_ggplotly$x$data[[i]]$name = gsub(',', ', ', obj_ggplotly$x$data[[i]]$name)
    }
  }
  return(obj_ggplotly)
}


#' Get figure title
#'
#' @param specimen specimen
#' @param display_type display unit type
#' @param du display unit
#'
#' @return Title
#'
#' @export
get_fig_title <- function(specimen, display_type, du){
  if(display_type != "Specimen"){
    #title = paste0("Specimen: ", specimen, " , ", display_type, ": ", du)
    title = paste0(specimen, " , ", du)
  }else{
    #title = paste0("Specimen: ", specimen)
    title = paste0(specimen)
  }
  return(title)
}


#' Concatenate file info
#'
#' @param v_i a vector of file numbers
#' @param v_name a vector of file names
#'
#' @return A vector of concatenated file info
#'
#' @export
cat_file_info <- function(i, name){
  return(paste(" ", sprintf("% 2d", i), name, sep="&emsp;"))
}

