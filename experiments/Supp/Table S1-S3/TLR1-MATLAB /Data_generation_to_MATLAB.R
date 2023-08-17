library(R.matlab)
source('./SimDataGeneration/true_signal_different.R')
source('./SimDataGeneration/data_generation_different.R')
n=1000

X_data_all <- array(0, c(64,64,n,50))
y_all_all <- array(0, c(5,n,50))
for(idx in 1:50){

set.seed(idx*100)

X_data=array(runif(64*64*n,0,1),c(64,64,n))
#X_data=array(rnorm(64*64*1000,0,1),c(64,64,n))
#X_data=array(rexp(64*64*1000,1),c(64,64,n))

X_datacross_reg_test=list()

X_datacross_reg_test[[1]]=X_data
X_datacross_reg_test[[2]]=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
X_datacross_reg_test[[3]]=X_data+adjust_nonlinear*0.5*cos(2*pi*X_data)

y_all=matrix(0,5,n)

y_all[1, ]=1+t(ctprod(BBprod[[1]],X_datacross_reg_test[[1]],2))+sigma_use[1]*rnorm(n,0,1)
y_all[2, ]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))+sigma_use[2]*rnorm(n,0,1)
y_all[3, ]=1+t(ctprod(BBprod[[3]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_test[[2]],2))+sigma_use[3]*rnorm(n,0,1)
y_all[4, ]=1+t(ctprod(BBprod[[5]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_test[[3]],2))+sigma_use[4]*rnorm(n,0,1)

y_all[5, ]=1+t(ctprod(BBprod_butterfly,X_datacross_reg_test[[2]],2))+sigma_use[5]*rnorm(n,0,1)


X_data_all[,,,idx]=X_data
y_all_all[,,idx]=y_all
}

writeMat("X_data_all_mat.mat", X_data_all=X_data_all)
writeMat("y_all_all_mat.mat", y_all_all=y_all_all)


