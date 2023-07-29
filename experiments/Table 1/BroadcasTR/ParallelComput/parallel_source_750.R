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
#source('./SimDataGeneration/true_signal_different.R')
#source('./SimDataGeneration/data_generation_different.R')
#source('./SimDataGeneration_new/20221018_get_inputs.R')
#source('./SimDataGeneration_new/20221018_get_BB.R')
#source('./SimDataGeneration_new/20221018_get_function.R')
#source('./SimDataGeneration_new/20221018_get_outputs.R', echo=FALSE)
#source('./SimDataGeneration_new/20221018_get_ise.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_inputs.R')
source('./SimDataGeneration_He/20230224_get_BB.R')
source('./SimDataGeneration_He/20230224_get_function.R')
source('./SimDataGeneration_He/20230224_get_outputs.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_ise.R', echo=FALSE)


R=c(1,2,3,4,5)
alpha=c(0,0.5,1)
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)


n_use=750


#rm(X_data_SNR)
#rm(X_datacross_reg_SNR)
#gc()

