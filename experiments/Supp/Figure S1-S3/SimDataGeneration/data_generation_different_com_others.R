
source('./SimDataGeneration/true_signal_different.R')

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
#X_data_SNR=array(runif(64*64*10000,0,1),c(64,64,10000))
#X_data_SNR=array(rnorm(64*64*10000,0,1),c(64,64,10000))
#X_data_SNR=array(rexp(64*64*10000,1),c(64,64,10000))
X_data_SNR = array(0, c(64,64,10000))
for(i in 1:10000){
    X_data_SNR[,,i]=as.matrix(runif(64)) %*% runif(64)
}



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




