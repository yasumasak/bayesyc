
library(ggplot2)
library(grid)
library(dplyr)

library(showtext) # To show Greek characters
showtext_auto() # To show Greek characters

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# Threshold for setting category
thre_alpha <- 5
thre_beta <- 0.1

prefix <- "figure_03."

v_sim_id <- c("lo", "mi", "hi")

v_color <- c(
  "peachpuff2", # 0.5
  "lightsalmon", # 1.0
  "hotpink2", # 2.5
  "red3", # 5.0
  "magenta4", # 7.5
  "purple4", # 10.0 
  "gray50"  # hierarchical
)

#------------------------------
# Plot normal distribution as a reference
v_x <- 10 ^ seq(-3.5, 3.5, 0.01)
v_sd <- c(0.5, 1, 2.5, 5, 7.5, 10)
v_p_all <- NULL
n_param <- length(v_sd)
for(sd in v_sd){
  v_p <- stats::dnorm(log(v_x), mean = 0, sd = sd)
  v_p_all <- c(v_p_all, v_p)
}
df_plot <- data.frame(x = rep(v_x, n_param), p = v_p_all,
                      sd = factor(rep(as.character(v_sd), each=length(v_x)),
                                  levels=as.character(v_sd)))

g <- ggplot(data=df_plot, aes(x=x, y=p, color=sd))
g <- g + theme_bw(base_size=10)
g <- g + scale_color_manual(values = v_color[1:6])
g <- g + geom_line(size=0.8, alpha=0.8)
g <- g + scale_x_log10(breaks=c(0.001, 0.01, 0.1, 1, 10, 100, 1000), labels=c("0.001", "0.01", "0.1", "1", "10", "100", "1000"))
g <- g + ylab("p(α)") + xlab("α")
g

# Output PDF
out_pdf_ <- paste0(prefix, "distribution_alpha.pdf")
pdf(file=out_pdf_, colormodel="cmyk", pointsize=7, paper="a4", width=0, height=0)
tryCatch({
  
  grid.newpage() # make new blank figure
  pushViewport(viewport(layout=grid.layout(4, 4)))
  vp.1  = viewport(layout.pos.row=1, layout.pos.col=1:2)
  print(g, vp=vp.1)
  
}, finally = { dev.off() }
)


#------------------------------
# Main analysis
for(data_type in c("count.0.01", "count.0.03")){
  
  E0 <- 10000
  base_dir <- paste0("../output/simulated_data/output.",data_type,".use_all")
  
  # models examined
  v_model <- c("bayesyc-i_nb.sigma_005",
               "bayesyc-i_nb.sigma_010",
               "bayesyc-i_nb.sigma_025",
               "bayesyc-i_nb.sigma_050",
               "bayesyc-i_nb.sigma_075",
               "bayesyc-i_nb.sigma_100",
               "bayesyc-i_nb.hierarchical"
              )
  
  # model names displayed
  v_model_name <- c("0.5",
                    "1.0",
                    "2.5",
                    "5.0",
                    "7.5",
                    "10.0",
                    "hierarchical"
                    )
  
  df_included_whole <- NULL
  df_param_whole <- NULL
  for(sim_id in v_sim_id){
    
    for(m in seq(v_model)){
      model <- v_model[m]
      model_name <- v_model_name[m]
      
      for(e1 in c(0.01, 0.33, 0.66) * E0){
      for(e2 in c(0.01, 0.33, 0.66) * E0){
      
      # From stats csv
      in_ <- file.path(base_dir, paste0("output_", model), 
                       paste0("stats.simulated_double.",data_type,".",sim_id,"_",e1,"_",e2,".",model,".csv")) # for using e1 and e2
      df_stats <- read.csv2(in_, header=T, sep=",", stringsAsFactors=F, check.names=F, quote="", comment.char="#")
      
      # model name
      df_stats$model <- model_name
      
      mx_param1 <- t(as.data.frame(strsplit(df_stats$Compound_1, "_\\w")))[, -1]
      colnames(mx_param1) <- c("C1", "E1", "h1", "Alpha12", "Beta")
      mx_param2 <- t(as.data.frame(strsplit(df_stats$Compound_2, "_\\w")))[, -1]
      colnames(mx_param2) <- c("C2", "E2", "h2", "Alpha21", "Beta")
      df_param <- as.data.frame(cbind(mx_param1, mx_param2[, c("C2", "E2", "h2", "Alpha21")]))
      rownames(df_param) <- seq(nrow(df_param))
      df_param <- as.data.frame(lapply(df_param, as.numeric))
      
      v_param_name <- c("Alpha12", "Alpha21", "Beta", "C1", "C2", "E1", "E2", "h1", "h2")
      df_included <- df_stats[, c("Compound_1", "Compound_2", "Success of fit", "model")]
      
      param_pattern <- paste(v_param_name, collapse="|")
      v_col_numeric <- grep(param_pattern, colnames(df_stats))
      df_stats[, v_col_numeric] <- lapply(df_stats[, v_col_numeric], as.numeric)
      
      #------------------------------
      # Check if the true value is included in 95% credible interval
      for(param in v_param_name){
        v_true_value <- df_param[, param]
        v_included <- ifelse(df_stats[, "Success of fit"] &
                             df_stats[, paste(param, "2.5%")] <= v_true_value & 
                             v_true_value <= df_stats[, paste(param, "97.5%")]
                             , 1, 0)
        df_included[, param] <- v_included
      }
      
      df_included$data_id <- sim_id
      df_included$data_type <- data_type
      df_included_whole <- rbind(df_included_whole, df_included)
      df_param_whole <- rbind(df_param_whole, df_param)
      
      }
      }
    }
    
  }
  
}
  

#--------------------------------------------------
# Output
for(data_type in c("all")){
#for(data_type in c("all", "count.0.01", "count.0.03")){
  v_true_alpha <- c(NA, 1, 10, 100)
  for(true_alpha in v_true_alpha){
    
    if(is.na(true_alpha)){
      out_pdf_ <- paste0(prefix, "sigma_alpha.estimates.95ci_included.simulated_double.",data_type,".bayesyc-i.pdf")
      
      # by data_type
      if(data_type == "all"){
        select_row <- rep(TRUE, nrow(df_param_whole))
      }else{
        select_row <- df_param_whole$data_type == data_type
      }
      
    }else{
      out_pdf_ <- paste0(prefix, "sigma_alpha.estimates.95ci_included.simulated_double.",data_type,".bayesyc-i.true_alpha_",true_alpha,".pdf")
      select_row <- df_param_whole$Alpha12 %in% c(true_alpha, 1/true_alpha) & 
                    df_param_whole$Alpha21 %in% c(true_alpha, 1/true_alpha)
      
      # by data_type
      if(data_type != "all"){
        select_row <- df_param_whole$Alpha12 %in% c(true_alpha, 1/true_alpha) & 
                      df_param_whole$Alpha21 %in% c(true_alpha, 1/true_alpha) & 
                      df_param_whole$data_type == data_type
      }
    }
    
    df_included_all <- df_included_whole[select_row, ]
    
    tb_examined <- table(df_included_all[ , c("model", "data_id")])
    n_examined <- unique(tb_examined)
    
    # Count included samples
    v_param_name_2 <- c("Alpha12", "Alpha21", "Beta", "C1", "C2", "E1", "E2")
    df_counts_included_all_2 <- df_included_all %>% select(model, data_id, Alpha12:E2) %>%
      group_by(model, data_id) %>%
      dplyr::summarise(across(Alpha12:E2, sum), .groups="keep")
    
    # Long format
    dfl_counts_included_all_2 <- df_counts_included_all_2 %>%
      tidyr::pivot_longer(cols=all_of(v_param_name_2), names_to="parameter", values_to="count")
    
    # Ratio
    dfl_counts_included_all_2$ratio <- dfl_counts_included_all_2$count / n_examined
    
    dfl_counts_included_all_2[dfl_counts_included_all_2$parameter=="Alpha12", "parameter"] <- "α12"
    dfl_counts_included_all_2[dfl_counts_included_all_2$parameter=="Alpha21", "parameter"] <- "α21"
    dfl_counts_included_all_2[dfl_counts_included_all_2$parameter=="Beta", "parameter"] <- "β"
    
    
    #--------------------------------------------------
    # Output PDF
    pdf(file=out_pdf_, colormodel="cmyk", pointsize=7, paper="a4", width=0, height=0)
    tryCatch({
    
      grid.newpage() # make new blank figure
      pushViewport(viewport(layout=grid.layout(5, 4)))
      
      #------------------------------
      # Plot for each simulation id
      for(i in seq(v_sim_id)){
        sim_id <- v_sim_id[i]
        df_plot <- dfl_counts_included_all_2[dfl_counts_included_all_2$data_id==sim_id, ]
        # Total
        v_mean_count <- tapply(df_plot$count, df_plot$model, mean)
        v_mean_ratio <- tapply(df_plot$ratio, df_plot$model, mean)
        df_avg <- data.frame(parameter="Total", model=names(v_mean_count), count=v_mean_count, data_id=sim_id, ratio=v_mean_ratio)
        v_param <- c(unique(df_plot$parameter), "Total")
        df_plot <- rbind(df_plot, df_avg)
        df_plot$model <- factor(df_plot$model, v_model_name)
        
        # Make plot
        g1 <- ggplot(df_plot, aes(x=factor(parameter, levels=v_param), y=as.numeric(ratio), fill=model))
        g1 <- g1 + theme_bw(base_size=8)
        g1 <- g1 + geom_bar(stat="identity", position=position_dodge(0.8), width = 0.7, alpha=1.0)
        g1 <- g1 + ggtitle(paste(sim_id))
        g1 <- g1 + scale_fill_manual(values=v_color)
        g1 <- g1 + theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="right")
        g1 <- g1 + coord_cartesian(ylim=c(0, 1))
        g1 <- g1 + xlab("Parameters") +ylab("Included in 95% CI")
        g1
        
        vp.1  = viewport(layout.pos.row=i, layout.pos.col=1:3)
        print(g1, vp=vp.1)
      }
      
      #------------------------------
      # Total across simulation ids
      dfl_counts_total <- dfl_counts_included_all_2 %>% group_by(parameter, model) %>%
        dplyr::summarise(total_count = sum(count), .groups = "keep")
      # Ratio
      dfl_counts_total$ratio <- dfl_counts_total$total_count / (n_examined * 3)
      
      # Total of total
      v_mean_count <- tapply(dfl_counts_total$total_count, dfl_counts_total$model, mean)
      v_mean_ratio <- tapply(dfl_counts_total$ratio, dfl_counts_total$model, mean)
      df_avg <- data.frame(parameter="Total", model=names(v_mean_count), total_count=v_mean_count, ratio=v_mean_ratio)
      df_plot_total <- rbind(dfl_counts_total, df_avg)
      df_plot_total$model <- factor(df_plot_total$model, v_model_name)
      
      # Make plot
      g1 <- ggplot(df_plot_total, aes(x=factor(parameter, levels=v_param), y=as.numeric(ratio), fill=model))
      g1 <- g1 + theme_bw(base_size=8)
      g1 <- g1 + geom_bar(stat="identity", position=position_dodge(0.8), width = 0.7, alpha=1.0)
      g1 <- g1 + ggtitle(paste("Total"))
      g1 <- g1 + scale_fill_manual(values=v_color)
      g1 <- g1 + theme(axis.text.x = element_text(angle=45, hjust=1), legend.position="right")
      g1 <- g1 + coord_cartesian(ylim=c(0, 1))
      g1 <- g1 + xlab("Parameters") +ylab("Included in 95% CI")
      g1
      
      vp.t  = viewport(layout.pos.row=i+1, layout.pos.col=1:3)
      print(g1, vp=vp.t)
      
    }, finally = { dev.off() }
    )
  
  }

}
