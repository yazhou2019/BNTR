#parallel computing for multiple replications


library(snowfall)
library(stringr)

sfInit(parallel=TRUE, cpus=10, type="SOCK")




wrapper <- function(idx) {
  # Output progress in worker logfile
  cat( "Current index: ", idx, "\n" )
  set.seed(idx)
  indexall = sample(1:774)
  BB_res = list()
  BB_idx=c()
  id_train = indexall[1:495]
  id_vali = indexall[496:618]
  id_test = indexall[619:774]
  res = validation_broadcasted_sparsetenreg(R,alpha,lambda,X_all[,,,id_train],y_all[id_train],X_all[,,,id_vali],y_all[id_vali],X_all[,,,id_test],y_all[id_test])
    BB_idx[1] = res$MSE_vali
    BB_idx[2]=res$MSE_pre
    
    BB_res[[1]] = BB_idx
    BB_res[[2]] = res$b_validation_test_lambda_R_nonlinear
    BB_res[[3]] = res$BB
    
  return(BB_res)
}

sfSource('./ParallelComput/parallel_source_ADNI.R')

sfClusterSetupRNG()

result <- sfLapply(1:10, wrapper)



sfStop()

setwd("./RealResults")
save(result,file = "ADNI40403.Rdata")
setwd("../")
