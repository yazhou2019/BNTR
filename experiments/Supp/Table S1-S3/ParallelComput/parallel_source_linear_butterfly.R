#parallel sourece
library(stringr)
library(rTensor)
library(glmnet)
library(Matrix)
library(MASS)
#parallel sourece
library(stringr)


source('./ComponentsLin/functions_needed.R')
source('./ComponentsLin/validation_result_linear.R')
source('./ComponentsLin/broadcasted_sparsetenreg_linear.R')

source('./SimDataGeneration/true_signal_different_butterfly.R')
source('./SimDataGeneration/data_generation_different_butterfly.R')


R=c(1,2,3,4,5)
alpha=c(0,0.5,1)
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)


n_use=250


rm(X_data_SNR)
rm(X_datacross_reg_SNR)
gc()


