#parallel computing for multiple replications


library(snowfall)
library(stringr)

sfInit(parallel=TRUE, cpus=5, type="SOCK")




wrapper <- function(idx) {
  # Output progress in worker logfile
  cat( "Current index: ", idx, "\n" )

  BB_idx=list()
  id_train = idtrain_matrix[idx,]
  id_vali = idtest_matrix[idx, 1:1000]
  id_test = idtest_matrix[idx, 1001:6000]
  
  res = validation_broadcasted_sparsetenreg(R,alpha,lambda,X_all[,,,id_train],y_all[id_train],X_all[,,,id_vali],y_all[id_vali],X_all[,,,id_test],y_all[id_test],num_knots=10)
    BB_idx[[1]] = res$MSE_vali
    BB_idx[[2]]=res$MSE_pre
    BB_idx[[3]]=res$b_validation_test_lambda_R_nonlinear
    index_best = which.min(res$b_validation_test_lambda_R_nonlinear[,2])
    BB_idx[[4]]=res$betabig[[index_best]]

  return(BB_idx)
}

sfSource('./ParallelComput/parallel_source.R')

sfClusterSetupRNG()

result <- sfLapply(1:10, wrapper)



sfStop()

setwd("./SimResults")
save(result,file = "monkey_1_new_K11_tuning_sim_new2.Rdata")
setwd("../")
