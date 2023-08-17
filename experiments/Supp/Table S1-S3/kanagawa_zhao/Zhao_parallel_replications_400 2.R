#parallel computing for multiple replications


library(snowfall)
sfInit(parallel=TRUE, cpus=10, type="SOCK")




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
  
 
  
  id_train=1:n_train
  id_vali=(1000-n_vali+1):1000
  id_test=851:1000 # test is useless here
  
  Xmat_list = list()
  Xmat_list[[1]] = array(0,c(64,64,n_train))
  Xmat_list[[2]] = array(0,c(64,64,n_train))
  library(rTensor)
  for(i in 1:n_train){
    tem = as.tensor(X_data[,,id_train[i]])
    Xmat_list[[1]][,,i]=k_unfold(tem,1)@data
    Xmat_list[[2]][,,i]=k_unfold(tem,2)@data
  }

  n_pre = 150
  set.seed(idx*1993)
  X_data_pre=array(runif(64*64*n_pre,0,1),c(64,64,n_pre))
  X_datacross_reg_test=list()
  X_datacross_reg_test[[1]]=X_data_pre
  X_datacross_reg_test[[2]]=X_data_pre+adjust_nonlinear*sin(2*pi*(X_data_pre-0.5)^2)
  X_datacross_reg_test[[3]]=X_data_pre+adjust_nonlinear*0.5*cos(2*pi*X_data_pre)
  y_all_pre=list()
  y_all_pre[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_test[[1]],2))+sigma_use[1]*rnorm(n_pre,0,1)
  y_all_pre[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))+sigma_use[2]*rnorm(n_pre,0,1)
  y_all_pre[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_test[[2]],2))+sigma_use[3]*rnorm(n_pre,0,1)
  y_all_pre[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_test[[3]],2))+sigma_use[4]*rnorm(n_pre,0,1)
  X_pre = X_data_pre
  
  Xmat_list_pre =list()
  Xmat_list_pre[[1]] = array(0,c(64,64,n_pre))
  Xmat_list_pre[[2]] = array(0,c(64,64,n_pre))
  for(i in 1:n_pre){
    tem = as.tensor(X_pre[,,i])
    Xmat_list_pre[[1]][,,i]=k_unfold(tem,1)@data
    Xmat_list_pre[[2]][,,i]=k_unfold(tem,2)@data
  }
  
  KXX_base_use = KXX_base_fun(Xmat_list, c(1,1,1))
  KXX_base_pre_use = KXX_base_fun_pre(Xmat_list_pre, Xmat_list, c(1,1,1))
  
  rm(X_datacross_reg_test)
  rm(X_pre)
  rm(X_data_pre)
  rm(X_data)
  ls()
  
  for(signal_i in 1:4){
    #res=validation_broadcasted_sparsetenreg(R,alpha,lambda,X_data[,,id_train],y_all[[signal_i]][id_train],X_data[,,id_vali],y_all[[signal_i]][id_vali],X_data[,,id_test],y_all[[signal_i]][id_test],num_knots=8)
    res_f = final_function(Xmat_list,Xmat_list_pre, y_all[[signal_i]][id_train], y_all_pre[[signal_i]], KXX_base_use, KXX_base_pre_use)
    BB_idx[[signal_i]] = res_f
    #BB_idx[[signal_i]]=res$BB
    #BB_idx[[signal_i+4]]=res$bb
    #BB_idx[[signal_i+8]]=res$MSE_pre
  }
  return(BB_idx)
}

sfSource('./ParallelComput/Zhao_source.R')

sfClusterSetupRNG()

result <- sfLapply(1:50, wrapper)


sfStop()
filename <- str_c("BNTR", n_use, "_Zhao400.Rdata")
setwd("./SimResults")
save(result,file = filename)
setwd("../")
