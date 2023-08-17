library(R.matlab)
library(stringr)
n_data = 500
distri = 'uniform'

BB_all = readMat(str_c("BB_all_",n_data,"_grid_butterfly.mat"))$BB.all
b0_all = readMat(str_c("b0_all_",n_data,"_grid_butterfly.mat"))$b0.all

result = list()
tem = list()
for(iter in 1:50){
  tem[[1]] = BB_all[,,5,iter]
  tem[[1+4]] = b0_all[5,iter]
  result[[iter]]=tem
}

savename = str_c("TLR1_", n_data, str_c("_linear_butterfly",".Rdata"))

save(file=savename, result)