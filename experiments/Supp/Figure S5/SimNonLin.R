#sample size
n_all = c(500, 750, 1000)
library(stringr)
# parallel computation, training the model using validation method
for(nn in c(1,3)){
    n=n_all[nn]
    filerep_name=str_c("./ParallelComput/parallel_replications_big_", n, ".R")
    source(filerep_name)
}

