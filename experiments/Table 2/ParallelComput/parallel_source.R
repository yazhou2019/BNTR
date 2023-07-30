#parallel sourece
#parallel sourece
library(stringr)
library(rTensor)
library(glmnet)
library(Matrix)
library(MASS)

source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')

R=c(1,2,3,4,5,6,7,8)
alpha=c(0,0.5,1)
#alpha=c(0)
lambda=c(0.001,0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1,2.5,5,7.5,10,25,50,75,100,250,500,750,1000)



# load the output

#load("/home/grad/zhou/monkey_1/y_output_1.Rdata")
#Y_big=Y_big[1:10000,1]
#rm(Y_output)


#set.seed(2018)
#id=sample(10000,10000)
#idtrain=id[1:4000]
#idtest=id[4001:10000]

#y_test=Y_big[idtest]
#y_train=Y_big[idtrain]

#X_test=X_big[,,,idtest]
#X_train=X_big[,,,idtrain]

#rm(Y_big)
#gc()


# load the input
load("/home/grad/zhou/monkey_1/X_test_1.Rdata")
load("/home/grad/zhou/monkey_1/X_train_1.Rdata")
# load("/home/grad/zhou/monkey_1/X_basis_test_monkey_1.Rdata")
# load("/home/grad/zhou/monkey_1/X_basis_train_monkey_1.Rdata")




# used for replication

#y_all=c()
#y_all[1:4000]=y_train
#y_all[4001:10000]=y_test
#X_all=array(0,c(64,10,10,6,10000))
X_all=array(0,c(64,10,10,10000))
 X_all[,,,1:4000]=X_train
 X_all[,,,4001:10000]=X_test
#X_all[,,,,1:4000]=X_basis_train
#X_all[,,,,4001:10000]=X_basis_test

#rm(X_basis_train)
#rm(X_basis_test)
rm(X_train)
rm(X_test)
gc()

#idtrain_matrix=matrix(0,20,4000)
#idtest_matrix=matrix(0,20,6000)
#for(iter in 1:20){
#    set.seed(100*iter+1018)
#    id=sample(10000,10000)
#    idtrain=id[1:4000]
#    idtest=id[4001:10000]

#    idtrain_matrix[iter,]=idtrain
#    idtest_matrix[iter,]=idtest
#}
load("/home/grad/zhou/gridsearch/monkey1/y_all.Rdata")
load("/home/grad/zhou/gridsearch/monkey1/idtest_matrix.Rdata")
load("/home/grad/zhou/gridsearch/monkey1/idtrain_matrix.Rdata")
