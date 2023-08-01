#In the theory, we only need the objective value is small enough (compared to the objective value for spline)




splinefunfit <- function(coeffit, knots, X){
  tildeX = array(0, c(64,64,7,1000))
  #tildex = cbind(1,x,x^2,x^3,(x- knots_all[[idx]][2])^3,(x-  knots_all[[idx]][3])^3,(x - knots_all[[idx]][4])^3)
  tildeX[,,1,] = 1
  tildeX[,,2,] = X
  tildeX[,,3,] = X^2
  tildeX[,,4,] = X^3
  tildeX[,,5,] = (X-knots[2])^3*((X- knots[2])^3>0)
  tildeX[,,6,] = (X-knots[3])^3*((X- knots[3])^3>0)
  tildeX[,,7,] = (X-knots[4])^3*((X- knots[4])^3>0)
  
  fX = array(0, c(64,64,1000))
  for(k in 1:7){
    fX = fX + tildeX[,,k,] * coeffit[k]
  }
  return(fX)
}


#parallel sourece
library(stringr)
library(rTensor)
library(glmnet)
library(Matrix)
library(MASS)

set.seed(2021)
n=500
method="BNTR"
kname = c("K3","K4","K5","K6","K7","K8","K9")
#kname = "K6"
#res_all_K = array(0, c(50,4,length(kname))
nk=4
path=str_c("./SimResultsKKK/",method,n, "_new_",kname[nk],".Rdata")
load(path)
  
  source('./ComponentsNonLin/functions_needed.R')
  source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
  source('./ComponentsNonLin/sequential_warmstart.R')
  source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')
  source('./SimDataGeneration/true_signal_different.R')
  source('./SimDataGeneration/data_generation_different.R')
  
  
  
  rm(X_data_SNR)
  rm(X_datacross_reg_SNR)
  gc()
  
  res_mat = matrix(0, 50, 4)
  res_mat_app = matrix(0, 50, 4)
  
  for(iter in 1:50){
    
  set.seed(iter* 100)

  n_use=500
  n_train=n_use * 0.8
  n_vali=n_use *0.2
  n = n_use
  
  id_train=1:n_train
  id_vali=(1000-n_vali+1):1000
  id_test=851:1000 # test 
  
  
  X_data=array(runif(64*64*1000,0,1),c(64,64,n))
  X_datacross_reg_test=list()
  X_datacross_reg_test[[1]]=X_data
  X_datacross_reg_test[[2]]=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
  X_datacross_reg_test[[3]]=X_data+adjust_nonlinear*0.5*cos(2*pi*X_data)
  
  y_all=list()
  
  y_all[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_test[[1]],2))+sigma_use[1]*rnorm(n,0,1)
  y_all[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))+sigma_use[2]*rnorm(n,0,1)
  y_all[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_test[[2]],2))+sigma_use[3]*rnorm(n,0,1)
  y_all[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_test[[3]],2))+sigma_use[4]*rnorm(n,0,1)
  
  #y_hat = list()
  
  knots = c(0, 0.25,0.5, 0.75,  1)
  x = 1:100000*0.00001
  #X_data+adjust_nonlinear
  signal_i = 1
  fftrue = list()
  fftrue[[1]] = x
  fftrue[[2]] = x + adjust_nonlinear *sin(2*pi*(x-0.5)^2)
  fftrue[[3]] = x + adjust_nonlinear*0.5*cos(2*pi*x)
  
  coef_truspline = list()
  
  #approxmate above functions
  
  tildex = cbind(1,x,x^2,x^3,(x- knots[2])^3*((x- knots[2])^3>0),(x-  knots[3])^3*((x-  knots[3])>0),(x - knots[4])^3*((x - knots[4])>0) )
  for(f_i in  1:3){
  coef_truspline[[f_i]] = solve(crossprod(tildex)) %*% t(tildex)%*%fftrue[[f_i]] 
  }
  

  yysplinepre = list()
  
  spline_reg = list()
  spline_reg[[1]] = splinefunfit(coeffit=coef_truspline[[1]], knots=knots, X=X_data)
  spline_reg[[2]] = splinefunfit(coeffit=coef_truspline[[2]], knots=knots, X=X_data)
  spline_reg[[3]] = splinefunfit(coeffit=coef_truspline[[3]], knots=knots, X=X_data)
  
  yysplinepre[[1]] = 1+t(ctprod(BBprod[[1]],spline_reg[[1]],2))
  yysplinepre[[2]] = 1+t(ctprod(BBprod[[2]],spline_reg[[2]],2))
  yysplinepre[[3]] = 1+t(ctprod(BBprod[[3]],spline_reg[[2]],2))+t(ctprod(BBprod[[4]],spline_reg[[2]],2))
  yysplinepre[[4]] = 1+t(ctprod(BBprod[[5]],spline_reg[[2]],2))+t(ctprod(BBprod[[6]],spline_reg[[3]],2)) 
  
  for(signal_i in 1:4){
    res_mat_app[iter,signal_i] = (sum((y_all[[signal_i]][id_train] - yysplinepre[[signal_i]][id_train])^2)/length(y_all[[signal_i]][ id_train]))
  }
  
  
  #num_knots <- nk+1
  #knots = quantile(c(X_data), probs = c(seq(0, 1, 1/(num_knots - 1))))
  knots = c(0, 0.25,0.5, 0.75,  1)
  tildePhiX_data = tildePhiX_trans(X_data[,,id_train], knots)
  
  for(signal_i in 1:4){
      BB=result[[iter]][[signal_i]]
      b0=result[[iter]][[signal_i+4]]
      y_hat= b0 + crossprod(matrix(tildePhiX_data, c(prod(dim(BB)), n_train)), as.vector(BB))
      res_mat[iter,signal_i] = sum((y_hat-y_all[[signal_i]][id_train])^2)/length(y_hat)
  }
  
  
  }
  
  