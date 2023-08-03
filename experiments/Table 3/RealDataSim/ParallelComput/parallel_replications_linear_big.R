#parallel computing for multiple replications


library(snowfall)
sfInit(parallel=TRUE, cpus=2, type="SOCK")



wrapper <- function(idx) {
    # Output progress in worker logfile
    cat( "Current index: ", idx, "\n" )
    
    BB_idx=list()
    id_train = idtrain_matrix[idx,]
    id_vali = idtest_matrix[idx, 1:1000]
    id_test = idtest_matrix[idx, 1001:6000]
    
    res = validation_result(R,alpha,lambda,X_all[,,,id_train],y_all[id_train],X_all[,,,id_vali],y_all[id_vali],X_all[,,,id_test],y_all[id_test])
    BB_idx[[1]] = res$MSE_vali
    BB_idx[[2]] = res$MSE_pre
    BB_idx[[3]] = res$BB
    BB_idx[[4]] = res$bb
    BB_idx[[5]] = res$knots
    
    return(BB_idx)
}

sfSource('./ParallelComput/parallel_source_linear.R')

sfClusterSetupRNG()

result <- sfLapply(1:10, wrapper)


sfStop()
filename <- str_c("TLR2_", "linear_new_sim_new2.Rdata")
setwd("./RealResults")
save(result,file = filename)
setwd("../")
