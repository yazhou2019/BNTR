
source('true_signal_different.R')

#noise level
sigma=0.1

adjust_nonlinear=0.3
#presentation
adjust_nonlinear=adjust_nonlinear*2








ctprod <- function(A,B,K){
  da <- dim(A)
  la <- length(da)
  db <- dim(B)
  lb <- length(db)
  Amat <- array(A, dim=c(prod(da[1:(la-K)]),prod(da[(la-K+1):la])))
  Bmat <- array(B, dim=c(prod(db[1:K]),prod(db[(K+1):lb])))
  Cmat <- Amat %*% Bmat
  C <- array(Cmat,dim=c(da[1:(la-K)],db[(K+1):lb]))
  return(C)
} 

set.seed(2018)
#data generation
#control SNR
X_data_SNR=array(runif(64*64*10000,0,1),c(64,64,10000))
#X_data_SNR=array(rnorm(64*64*10000,0,1),c(64,64,10000))
#X_data_SNR=array(rexp(64*64*10000,1),c(64,64,10000))

X_datacross_reg_SNR=list()
X_datacross_reg_SNR[[1]]=X_data_SNR
X_datacross_reg_SNR[[2]]=X_data_SNR+adjust_nonlinear*sin(2*pi*(X_data_SNR-0.5)^2)
X_datacross_reg_SNR[[3]]=X_data_SNR+adjust_nonlinear*0.5*cos(2*pi*X_data_SNR)


BBprod=list()
SNR_index=c()
for(signal_i in 1:6){
    BBprod[[signal_i]]=array(0,c(1,dim(BB_signal_all[[signal_i]])))
    BBprod[[signal_i]][1,,]=BB_signal_all[[signal_i]]
}


#BB_signal ++++++++++++++all
y_SNR=list()

y_SNR[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_SNR[[1]],2))
y_SNR[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_SNR[[2]],2))
y_SNR[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_SNR[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_SNR[[2]],2))
y_SNR[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_SNR[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_SNR[[3]],2))


sigma_use=c()
for(signal_i in 1:4){
  sigma_use[signal_i]=sigma*sqrt(var(y_SNR[[signal_i]]))
}



#test data set
#set.seed(77840)
#n_test=300

#X_data_test=array(runif(64*64*1000,0,1),c(64,64,n_test))
#X_data_test=array(rnorm(64*64*1000,0,1),c(64,64,n_test))
#X_data_test=array(rexp(64*64*1000,1),c(64,64,n_test))

#X_datacross_reg_test=list()

#X_datacross_reg_test[[1]]=X_data_test
#X_datacross_reg_test[[2]]=X_data_test+adjust_nonlinear*sin(2*pi*(X_data_test-0.5)^2)
#X_datacross_reg_test[[3]]=X_data_test+adjust_nonlinear*0.5*cos(2*pi*X_data_test)

#y_test=list()

#y_test[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_test[[1]],2))+sigma_use[1]*rnorm(n_test,0,1)
#y_test[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))+sigma_use[2]*rnorm(n_test,0,1)
#y_test[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_test[[2]],2))+sigma_use[3]*rnorm(n_test,0,1)
#y_test[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_test[[3]],2))+sigma_use[4]*rnorm(n_test,0,1)



#train data set
#set.seed(2018)
#n_size_sequnce=c(700)

#n=n_size_sequnce[1]
#X_data_train=array(runif(64*64*n,0,1),c(64,64,n))
#X_data_train=array(rnorm(64*64*n,0,1),c(64,64,n))
#X_data_train=array(rexp(64*64*n,1),c(64,64,n))

#X_datacross_reg_train=list()
#X_datacross_reg_train[[1]]=X_data_train
#X_datacross_reg_train[[2]]=X_data_train+adjust_nonlinear*sin(2*pi*(X_data_train-0.5)^2)
#X_datacross_reg_train[[3]]=X_data_train+adjust_nonlinear*0.5*cos(2*pi*X_data_train)

#y_train=list()

#y_train[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_train[[1]],2))+sigma_use[1]*rnorm(n,0,1)
#y_train[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_train[[2]],2))+sigma_use[2]*rnorm(n,0,1)
#y_train[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_train[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_train[[2]],2))+sigma_use[3]*rnorm(n,0,1)
#y_train[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_train[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_train[[3]],2))+sigma_use[4]*rnorm(n,0,1)


#rm(X_datacross_reg_train)
#rm(X_datacross_reg_SNR)
#rm(X_datacross_reg_test)
#rm(X_data_SNR)
#gc()

