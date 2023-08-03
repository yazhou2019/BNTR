
#parallel computing for multiple replications


library(snowfall)
sfInit(parallel=TRUE, cpus=10, type="SOCK")


#library(glmnet)
#source('./ENetComponents/ENet_source.R')



# # load the input
# load("./data/X_test_1.Rdata")
# load("./data/X_train_1.Rdata")
# # used for replication
# X_all=array(0,c(64,10,10,10000))
# X_all[,,,1:4000]=X_train
# X_all[,,,4001:10000]=X_test
# rm(X_train)
# rm(X_test)
# gc()
# load("./data/y_all_sim_new.Rdata")
# load("./data/idtest_matrix.Rdata")
# load("./data/idtrain_matrix.Rdata")




wrapper <- function(idx) {

BB_idx=list()
result=list()
#all the data set
set.seed(idx*10)
BB_idx=c()
id_train = idtrain_matrix[idx,]
id_vali = idtest_matrix[idx, 1:1000]
id_test = idtest_matrix[idx, 1001:6000]
res=validation_ENet(alpha,lambda,X_all[,,,id_train],y_all[id_train],X_all[,,,id_vali],y_all[id_vali],X_all[,,,id_test],y_all[id_test])
 

BB_idx[[1]]=matrix(res$BB, c(64,10,10))
BB_idx[[1+5]]=res$bb
BB_idx[[1+10]]=res$MSE_pre
BB_idx[[1+15]] = "ENet"
  
result[[idx]]=BB_idx

return(result)
}

sfSource('./ENetComponents/ENet_source.R')

sfClusterSetupRNG()

result <- sfLapply(1:10, wrapper)


filename <- str_c("ENet_real_20230103_2.Rdata")
setwd("./RealResults")
save(result,file = filename)
setwd("../")


