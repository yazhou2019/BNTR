#sample size
n_all = c(500, 750, 1000)
library(stringr)
# parallel computation, training the model using validation method
for(nn in 1){
    n=n_all[nn]
    filerep_name=str_c("./ParallelComput/parallel_replications_Zhao_", n, "_new.R")
    source(filerep_name)
}


# calculate the ISE
#for (nn in 1:3){
#n=n_all[nn]
#source('./CollectISEs/results_2_ISE_nonlinear.R')
#}
