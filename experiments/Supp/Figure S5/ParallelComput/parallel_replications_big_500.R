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
  
  X_data=array(runif(64*64*1000,0,1),c(64,64,n))
  #X_data=array(rnorm(64*64*1000,0,1),c(64,64,n))
  #X_data=array(rexp(64*64*1000,1),c(64,64,n))
  
  X_datacross_reg_test=list()
  
  X_datacross_reg_test[[1]]=X_data
  X_datacross_reg_test[[2]]=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
  X_datacross_reg_test[[3]]=X_data+adjust_nonlinear*0.5*cos(2*pi*X_data)
  
  y_all=list()
  
  y_all[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_test[[1]],2))+sigma_use[1]*rnorm(n,0,1)
  y_all[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))+sigma_use[2]*rnorm(n,0,1)
  y_all[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_test[[2]],2))+sigma_use[3]*rnorm(n,0,1)
  y_all[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_test[[3]],2))+sigma_use[4]*rnorm(n,0,1)
  
  
  
  
  for(signal_i in 1:4){
    
    id_train=1:n_train
    id_vali=(1000-n_vali+1):1000
    id_test=851:1000 # test is useless here
    
    res=validation_broadcasted_sparsetenreg(R,alpha,lambda,X_data[,,id_train],y_all[[signal_i]][id_train],X_data[,,id_vali],y_all[[signal_i]][id_vali],X_data[,,id_test],y_all[[signal_i]][id_test])
    BB_idx[[signal_i]]=res$BB
    BB_idx[[signal_i+4]]=res$bb
    BB_idx[[signal_i+8]]=res$MSE_pre
    BB_idx[[signal_i+12]]=res$b_validation_test_lambda_R_nonlinear
    index_best = which.min(res$b_validation_test_lambda_R_nonlinear[,2])
    BB_idx[[signal_i+16]]=res$betabig[[index_best]]
  }
  return(BB_idx)
}

sfSource('./ParallelComput/parallel_source_500.R')

sfClusterSetupRNG()

result <- sfLapply(1:50, wrapper)


sfStop()
filename <- str_c("BNTR", n_use, "_new_obj.Rdata")
setwd("./SimResults")
save(result,file = filename)
setwd("../")
