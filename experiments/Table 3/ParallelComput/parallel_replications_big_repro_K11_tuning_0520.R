#parallel computing for multiple replications


library(snowfall)
library(stringr)

sfInit(parallel=TRUE, cpus=5, type="SOCK")




wrapper <- function(idx) {
  # Output progress in worker logfile
  cat( "Current index: ", idx, "\n" )
  
  BB_idx=list()
  #data0520 = list(input_mat=mat_res$input_mat, ouput_mat=mat_res$output_mat, y_wo_error=y_wo_error, y_error=y_error)
  input_mat = data0520$input_mat
  output_mat =  data0520$ouput_mat
  y_wo_error = data0520$y_wo_error
  y_error = data0520$y_error
  
  id_input = input_mat[idx,]
  id_output = output_mat[idx,]
  
  y_all = y_wo_error[id_input] + y_error[id_output] 

  
  id_train = id_input[1:4000]
  id_vali = id_input[4001:5000]
  id_test = id_input[5001:10000]
  
  res = validation_broadcasted_sparsetenreg(R,alpha,lambda,X_all[,,,id_train],y_all[id_train],X_all[,,,id_vali],y_all[id_vali],X_all[,,,id_test],y_all[id_test],num_knots=10)
  BB_idx[[1]] = res$MSE_vali
  BB_idx[[2]]=res$MSE_pre
  BB_idx[[3]]=res$b_validation_test_lambda_R_nonlinear
  index_best = which.min(res$b_validation_test_lambda_R_nonlinear[,2])
  BB_idx[[4]]=res$betabig[[index_best]]
  
  return(BB_idx)
}

sfSource('./ParallelComput/parallel_source_0520.R')

sfClusterSetupRNG()

result <- sfLapply(1:10, wrapper)



sfStop()

setwd("./RealResults")
save(result,file = "monkey_1_new_K11_tuning_sim_new_0520.Rdata")
setwd("../")
