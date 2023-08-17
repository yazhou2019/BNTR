library(stringr)
source('./ComponentsLin/functions_needed.R')
source('./ComponentsLin/validation_result_linear.R')
source('./ComponentsLin/broadcasted_sparsetenreg_linear.R')

source('./SimDataGeneration/true_signal_different_butterfly.R')
source('./SimDataGeneration/data_generation_different_butterfly.R')




source('./ENet_butterfly/Validation_Enet.R')

alpha=c(0,0.5,1)
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)


BB_idx=list()
result=list()

for(idx in 1:50){


#all the data set
set.seed(idx*100)
n=1000

n_use=1000
n_train=n_use*0.8
n_vali=n_use*0.2

X_data=array(runif(64*64*1000,0,1),c(64,64,n))
#X_data=array(rnorm(64*64*1000,0,1),c(64,64,n))
#X_data=array(rexp(64*64*1000,1),c(64,64,n))




X_datacross_reg_test = X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)

y_all = 1 + t(ctprod(BBprod,X_datacross_reg_test,2))+sigma_use*rnorm(n,0,1)

for(signal_i in 1){
  
  id_train=1:n_train
  id_vali=(n-n_vali+1):n
  id_test=851:1000
  
  res=validation_ENet(alpha,lambda,X_data[,,id_train],y_all[id_train],X_data[,,id_vali],y_all[id_vali],X_data[,,id_test],y_all[id_test])
  BB_idx[[signal_i]]=matrix(res$BB, 64,64)
  BB_idx[[signal_i+4]]=res$bb
  BB_idx[[signal_i+8]]=res$MSE_pre

}


result[[idx]]=BB_idx

}

n_data = n_use

savename = str_c("./SimResults/ENet_", n_data, str_c("_linear_butterfly",".Rdata"))

save(file=savename, result)
