library(stringr)
library(stringi)
library(glmnet)
source('./ENetComponents/Validation_Enet.R')
#source('./SimDataGeneration_new/20221018_get_ise.R')

alpha=c(0,0.5,1)
#alpha=c(0)
lambda=c(0.001,0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1,2.5,5,7.5,10,25,50,75,100,250,500,750,1000)



# load the input
load("./data/X_test_1.Rdata")
load("./data/X_train_1.Rdata")
# used for replication
X_all=array(0,c(64,10,10,10000))
X_all[,,,1:4000]=X_train
X_all[,,,4001:10000]=X_test
rm(X_train)
rm(X_test)
gc()
load("./data/y_all_sim_new.Rdata")
load("./data/idtest_matrix.Rdata")
load("./data/idtrain_matrix.Rdata")
