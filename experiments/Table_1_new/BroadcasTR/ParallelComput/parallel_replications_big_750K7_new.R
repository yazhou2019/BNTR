#parallel computing for multiple replications


library(snowfall)
sfInit(parallel=TRUE, cpus=17, type="SOCK")




wrapper <- function(idx) {
  # Output progress in worker logfile
  cat( "Current index: ", idx, "\n" )
  BB_idx=list()
  
  #all the data set
  set.seed(idx*100)
  n=1000
  
  ##########################user define in the parallel_source########################################
  n_use=n_use
  ##########################user define in the parallel_source########################################
  n_train=n_use * 0.8
  n_vali=n_use *0.2
  
  #X_data=array(runif(64*64*1000,0,1),c(64,64,n))
  #X_data=array(rnorm(64*64*1000,0,1),c(64,64,n))
  #X_data=array(rexp(64*64*1000,1),c(64,64,n))
  #X_datacross_reg_test=list()
  #X_datacross_reg_test[[1]]=X_data
  #X_datacross_reg_test[[2]]=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
  #X_datacross_reg_test[[3]]=X_data+adjust_nonlinear*0.5*cos(2*pi*X_data)
  
  n_all = n_train + n_vali
  
  load("./SimResults/Simul_n1000_rep50_final_fix_new_editor.Rdata")
  
  start_time <- Sys.time()
  
  
  BB_list = data_all[[51]]
  v_list = data_all[[52]]
  data_all = data_all[[idx]]
  X_data = data_all$X
  y_all = data_all$y
  Max_iter = 0
  for(signal_i in 1:5){
    
    id_train=1:n_train
    id_vali=(1000-n_vali+1):1000
    id_test=851:1000 # test is useless here
    
    res=validation_broadcasted_sparsetenreg(R,alpha,lambda,X_data[,,id_train],y_all[[signal_i]][id_train],X_data[,,id_vali],y_all[[signal_i]][id_vali],X_data[,,id_test],y_all[[signal_i]][id_test],num_knots=6)
    BB_idx[[signal_i]]=res$BB
    BB_idx[[signal_i+5]]=res$bb
    BB_idx[[signal_i+10]]=res$MSE_pre
    BB_idx[[signal_i+15]] = res$knots_used
    BB_idx[[signal_i+20]] = BB_list[[signal_i]]
    BB_idx[[signal_i+25]] = v_list[[signal_i]]
    Max_iter = max(Max_iter, res$max_iter)
  }
  
  
  end_time <- Sys.time()
  max_mem <- max(mem_used() - baseline_mem)
  max_mem_in_mb <- max_mem / (1024 * 1024)
  
  BB_idx[[31]] = end_time - start_time
  BB_idx[[32]] = max_mem_in_mb
  BB_idx[[33]] = Max_iter
  
  return(BB_idx)
}

sfSource('./ParallelComput/parallel_source_750.R')

sfClusterSetupRNG()

result <- sfLapply(1:50, wrapper)


sfStop()
filename <- str_c("BNTR_", n_use, "_20230918.Rdata")
setwd("./SimResults")
save(result,file = filename)
setwd("../")
