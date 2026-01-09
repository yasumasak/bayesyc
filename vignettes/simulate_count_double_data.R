

pow <- function(x, y) x^y
  
musyc <- function(d1, d2, E0, E1, E2, E3, C1, C2, h1, h2, a12, a21, g12, g21, r1r, r2r){
    
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

opts <- list()

opts$v_conc <- c(0, 0.01, 0.1, 1, 10, 100)

opts$v_ic50 <- c(1, 0.05, 20)
opts$v_sim_id <- c("mi", "lo", "hi")

opts$v_alpha <- c(0.01, 0.1, 1, 10, 100)

opts$v_slope <- c(1)
opts$v_beta <- c(-0.4, 0, 0.4, 0.8, 1.6)

opts$E0 <- 10000
#opts$s_y <- 0.03 # CV=0.173
opts$s_y <- 0.01 # CV=0.1
opts$n_rep <- 3
opts$g12 <- 1
opts$g21 <- 1
opts$r1r <- 0.0001
opts$r2r <- 0.0001

opts$v_min <- c(0.01, 0.33, 0.66) * opts$E0
dir_out <- "../input/simulated_data"

dir.create(dir_out, showWarnings=FALSE, mode="0755")

set.seed(777)

for (i in seq(opts$v_sim_id)){
  sim_id <- opts$v_sim_id[i]
  C1 <- opts$v_ic50[i]
  C2 <- opts$v_ic50[i]
  count <- 0
  
  for (h1 in opts$v_slope){
  for (h2 in opts$v_slope){
    for (E1 in opts$v_min){
    for (E2 in opts$v_min){
      
      # initialize output vectors
      v_out_comp1 <- c()
      v_out_comp2 <- c()
      v_out_comb <- c()
      v_out_d1 <- c()
      v_out_d2 <- c()
      v_out_C1 <- c()
      v_out_C2 <- c()
      v_out_E1 <- c()
      v_out_E2 <- c()
      v_out_h1 <- c()
      v_out_h2 <- c()
      v_out_alpha12 <- c()
      v_out_alpha21 <- c()
      v_out_beta <- c()
      v_out_y <- c()
      
      for (a12 in opts$v_alpha){
      for (a21 in opts$v_alpha){
        for (beta in opts$v_beta){
          if(min(E1, E2) - beta * (opts$E0 - min(E1, E2)) < 0) next
          if(min(E1, E2) - beta * (opts$E0 - min(E1, E2)) > 0.67 * opts$E0) next
          E3 <- min(E1, E2) - beta * (opts$E0 - min(E1, E2))
          # sampling at each drug concentration
          for (d1 in opts$v_conc){
          for (d2 in opts$v_conc){
            mu <- musyc(d1, d2, opts$E0, E1, E2, E3, C1, C2, h1, h2, a12, a21, opts$g12, opts$g21, opts$r1r, opts$r2r)
            sim_y <- rnbinom(n=opts$n_rep, size=1/opts$s_y, mu=mu) # size = 1 / alpha
            v_out_comp1 <- c(v_out_comp1, rep(paste0("comp1_C",C1,"_E",E1,"_h",h1,"_A",a12,"_B",beta), opts$n_rep))
            v_out_comp2 <- c(v_out_comp2, rep(paste0("comp2_C",C2,"_E",E2,"_h",h2,"_A",a21,"_B",beta), opts$n_rep))
            v_out_comb <- c(v_out_comb, rep(paste0("C",C1,"_E",E1,"_h",h1,"_C",C2,"_E",E2,"_h",h2,"_A",a12,"_A",a21,"_B",beta), opts$n_rep))
            v_out_d1 <- c(v_out_d1, rep(d1, opts$n_rep))
            v_out_d2 <- c(v_out_d2, rep(d2, opts$n_rep))
            v_out_C1 <- c(v_out_C1, rep(C1, opts$n_rep))
            v_out_C2 <- c(v_out_C2, rep(C2, opts$n_rep))
            v_out_E1 <- c(v_out_E1, rep(E1, opts$n_rep))
            v_out_E2 <- c(v_out_E2, rep(E2, opts$n_rep))
            v_out_h1 <- c(v_out_h1, rep(h1, opts$n_rep))
            v_out_h2 <- c(v_out_h2, rep(h2, opts$n_rep))
            v_out_alpha12 <- c(v_out_alpha12, rep(a12, opts$n_rep))
            v_out_alpha21 <- c(v_out_alpha21, rep(a21, opts$n_rep))
            v_out_beta <- c(v_out_beta, rep(beta, opts$n_rep))
            v_out_y <- c(v_out_y, sim_y)
          }
          }
          count <- count + 1
        }
      }
      }
      
      df_sim <- data.frame(compound_1=v_out_comp1, compound_2=v_out_comp2, combination_name=v_out_comb,
                           conc_1=v_out_d1, conc_2=v_out_d2, y=v_out_y,
                           alpha12=v_out_alpha12, alpha21=v_out_alpha21, beta=v_out_beta,
                           C1=v_out_C1, C2=v_out_C2, E1=v_out_E1, E2=v_out_E2, h1=v_out_h1, h2=v_out_h2) 
      write.table(df_sim, file=paste0(dir_out,"/simulated_count_data.double.",opts$s_y,".",sim_id,"_",E1,"_",E2,".tsv"), sep="\t", row.names=F, col.names=T)

      # Print simulation parameters
      out_opts_ <- paste0(dir_out,"/simulated_count_data.double.",opts$s_y,".",sim_id,"_",E1,"_",E2,".opts.txt")
      sink(out_opts_)
      for(p in names(opts)){
        cat(p, "\t", paste(opts[[p]], collapse=","), "\n", sep="")
      }
      sink()
      
    }
    }
  }
  }
  print(paste(sim_id, count))
}

