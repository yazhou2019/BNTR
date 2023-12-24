#parallel sourece
library(stringr)


library(pryr)

baseline_mem <- mem_used()

source('./ComponentsLin/functions_needed.R')
source('./ComponentsLin/validation_result_linear.R')
source('./ComponentsLin/broadcasted_sparsetenreg_linear.R')
#source('./SimDataGeneration/true_signal_different.R')
#source('./SimDataGeneration/data_generation_different.R')

# source('./SimDataGeneration_new/20221018_get_inputs.R')
# source('./SimDataGeneration_new/20221018_get_BB.R')
# source('./SimDataGeneration_new/20221018_get_function.R')
# source('./SimDataGeneration_new/20221018_get_outputs.R', echo=FALSE)
# source('./SimDataGeneration_new/20221018_get_ise.R', echo=FALSE)

source('./SimDataGeneration_He/20230224_get_inputs.R')
source('./SimDataGeneration_He/20230224_get_BB.R')
source('./SimDataGeneration_He/20230224_get_function.R')
source('./SimDataGeneration_He/20230224_get_outputs.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_ise.R', echo=FALSE)


#####################################user define in the parallel_source_linear ######################
R=c(1,2,3,4,5)
#R=c(1,2,3,4,5,6,7,8)
alpha=c(0,0.5,1)
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)

n_use=1000
#####################################user define in the parallel_source_linear ######################




#set.seed(2019)
#id_matrix=matrix(0,100,1000)
#for(iter in 1:100){
#id_matrix[iter,]=sample(1000,1000)
#}




#signal_i=2
#nsequence=c(250,1000)
#num_n=length(nsequence)

