#parallel sourece
library(stringr)


source('./ComponentsLin/functions_needed.R')
source('./ComponentsLin/validation_result_linear.R')
source('./ComponentsLin/broadcasted_sparsetenreg_linear.R')
source('./SimDataGeneration/true_signal_different.R')
source('./SimDataGeneration/data_generation_different.R')

#####################################user define in the parallel_source_linear ######################
R=c(1,2,3,4,5)
# R=1
alpha=c(0,0.5,1)
#  alpha=1
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)
# lambda=c(0.01)
n_use=250
#####################################user define in the parallel_source_linear ######################


rm(X_data_SNR)
rm(X_datacross_reg_SNR)
gc()


#set.seed(2019)
#id_matrix=matrix(0,100,1000)
#for(iter in 1:100){
#id_matrix[iter,]=sample(1000,1000)
#}




#signal_i=2
#nsequence=c(250,1000)
#num_n=length(nsequence)

