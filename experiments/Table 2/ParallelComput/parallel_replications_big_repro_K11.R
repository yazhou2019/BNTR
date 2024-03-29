#parallel computing for multiple replications


library(snowfall)
library(stringr)

sfInit(parallel=TRUE, cpus=5, type="SOCK")




wrapper <- function(idx) {
  # Output progress in worker logfile
  cat( "Current index: ", idx, "\n" )

  BB_idx=c()
  id_train = idtrain_matrix[idx,]
  id_vali = idtest_matrix[idx, 1:1000]
  id_test = idtest_matrix[idx, 1001:6000]
  
  res = validation_broadcasted_sparsetenreg(R,alpha,lambda,X_all[,,,id_train],y_all[id_train],X_all[,,,id_vali],y_all[id_vali],X_all[,,,id_test],y_all[id_test],num_knots=10)
    BB_idx[1] = res$MSE_vali
    BB_idx[2]=res$MSE_pre

  return(BB_idx)
}

sfSource('./ParallelComput/parallel_source_repro.R')

sfClusterSetupRNG()

result <- sfLapply(1:10, wrapper)



sfStop()

setwd("./RealResults")
save(result,file = "monkey_1_new_K11.Rdata")
setwd("../")
