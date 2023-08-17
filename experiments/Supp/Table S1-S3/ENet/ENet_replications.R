source('~/Desktop/Validation_Enet.R')

source('~/Desktop/simulation for paper/Simulation_paper_new/ENet/data_generation_different.R')
alpha=c(0,0.5,1)
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)


BB_idx=list()
result=list()

for(idx in 1:50){


#all the data set
set.seed(idx*100)
n=1000

n_use=500
n_train=n_use*0.8
n_vali=n_use*0.2

X_data=array(runif(64*64*1000,0,1),c(64,64,n))
#X_data=array(rnorm(64*64*1000,0,1),c(64,64,n))
#X_data=array(rexp(64*64*1000,1),c(64,64,n))

X_datacross_reg_test=list()

X_datacross_reg_test[[1]]=X_data
X_datacross_reg_test[[2]]=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
X_datacross_reg_test[[3]]=X_data+adjust_nonlinear*0.5*cos(2*pi*X_data)

y_all=list()

y_all[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_test[[1]],2))+sigma_use[1]*rnorm(n,0,1)
y_all[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))+sigma_use[2]*rnorm(n,0,1)
y_all[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_test[[2]],2))+sigma_use[3]*rnorm(n,0,1)
y_all[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_test[[3]],2))+sigma_use[4]*rnorm(n,0,1)


for(signal_i in 1:4){
  
  id_train=1:n_train
  id_vali=(n-n_vali+1):n
  id_test=851:1000
  
  res=validation_ENet(alpha,lambda,X_data[,,id_train],y_all[[signal_i]][id_train],X_data[,,id_vali],y_all[[signal_i]][id_vali],X_data[,,id_test],y_all[[signal_i]][id_test])
  BB_idx[[signal_i]]=matrix(res$BB, 64,64)
  BB_idx[[signal_i+4]]=res$bb
  BB_idx[[signal_i+8]]=res$MSE_pre

}


result[[idx]]=BB_idx

}