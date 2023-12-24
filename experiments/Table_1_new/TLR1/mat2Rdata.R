#setwd("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper")
setwd("~/Desktop/Research/JRSSB/换editor上传/Table_1")
library(R.matlab)
library(stringr)
library(stringi)

n_use = 500

BBname = str_c("./TLR1/BB_all_", n_use, "_grid.mat")
b0name = str_c("./TLR1/b0_all_", n_use, "_grid.mat")
MSEname = str_c("./TLR1/MSE_all_", n_use, "_grid.csv")
  
BB_all = readMat(BBname)$BB.all
b0_all = readMat(b0name)$b0.all
MSE_all = as.matrix(read.csv(MSEname, header = FALSE))


print(dim(BB_all))
print(dim(b0_all))

result = list()
for(iter in 1:50){
  result[[iter]] = list()
  for(signal_i in 1:5){
    result[[iter]][[signal_i]] = BB_all[,,signal_i, iter]
    result[[iter]][[signal_i + 5]] = b0_all[signal_i, iter]
    result[[iter]][[signal_i + 10]] = MSE_all[iter, signal_i]
  }
  
}

filename <- str_c("./TLR1/TLR1_", n_use,"_20230918.Rdata")
save(result,file = filename)


# BB = result[[idx]][[signal_i]]
# bb = result[[idx]][[signal_i+5]]
# MSE_pre = result[[idx]][[signal_i+10]]
# knots_used = result[[idx]][[signal_i+15]]
# #BBv = result[[idx]][[signal_i+20]]
# vv = result[[idx]][[signal_i+25]]
# ise_mat[1, signal_i] = final_ise_linear(b0=bb, BB=BB, BBv=BB_list, v =vv, case=signal_i, precision=1000, knots = knots_used, order = 4)$ise
# ise_mat[2, signal_i] = MSE_pre
