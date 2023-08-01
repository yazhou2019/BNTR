#parallel sourece
library(stringr)
library(rTensor)
library(glmnet)
library(Matrix)
library(MASS)


#In the theory, we only need the objective value is small enough (compared to the objective value for spline)


svd_to_cp <- function(BB,R=4){
  tem = svd(BB)
  bbeta = list()
  bbeta[[1]] = tem$u %*% diag(tem$d)
  bbeta[[2]] = tem$v
  
  return(bbeta)
}


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


set.seed(2021)
source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')
source('./SimDataGeneration/true_signal_different.R')
source('./SimDataGeneration/data_generation_different.R')
rm(X_data_SNR)
rm(X_datacross_reg_SNR)
gc()


n=1000
n_use = 500
path = str_c("./SimResults/BNTR",n_use, "_new_tuning_par.Rdata")
load(path)


res_mat = matrix(0, 50, 4)
res_mat_app = matrix(0, 50, 4)

for(iter in 1:50){
  
  set.seed(iter* 100)
  
  n_use=n_use
  n_train=n_use * 0.8
  n_vali=n_use *0.2
  n = n_use
  
  id_train=1:n_train
  id_vali=(1000-n_vali+1):1000
  id_test=851:1000 # test 
  
  
  X_data=array(runif(64*64*1000,0,1),c(64,64,n))
  X_datacross_reg_test=list()
  X_datacross_reg_test[[2]]=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
  
  y_all=list()
  y_all[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))+sigma_use[2]*rnorm(n,0,1)
  
  
  
  #knots = c(0, 0.25,0.5, 0.75,  1)
  num_knots = 5
  knots = quantile(c(X_data[,,id_train]), probs = c(seq(0, 1, 1/(num_knots - 1))))
  
  x = 1:100000*0.00001
  #X_data+adjust_nonlinear
  fftrue = list()
  fftrue[[2]] = x + adjust_nonlinear *sin(2*pi*(x-0.5)^2)
  
  coef_truspline = list()
  
  #approxmate above functions
  tildex = cbind(1,x,x^2,x^3,(x- knots[2])^3*((x- knots[2])^3>0),(x-  knots[3])^3*((x-  knots[3])>0),(x - knots[4])^3*((x - knots[4])>0) )
  for(f_i in  2){
    coef_truspline[[f_i]] = solve(crossprod(tildex)) %*% t(tildex)%*%fftrue[[f_i]] 
  }
  
  
  yysplinepre = list()
  
  spline_reg = list()
  spline_reg[[2]] = splinefunfit(coeffit=coef_truspline[[2]], knots=knots, X=X_data)
  yysplinepre[[2]] = 1+t(ctprod(BBprod[[2]],spline_reg[[2]],2))
  
  #for(signal_i in 2){
  #  res_mat_app[iter,signal_i] = (sum((y_all[[signal_i]][id_train] - yysplinepre[[signal_i]][id_train])^2)/length(y_all[[signal_i]][ id_train]))
  # }
  
  
  
  tildePhiX_data = tildePhiX_trans(X_data[,,id_train], knots)
  
  
  for(signal_i in 2){
    BB=result[[iter]][[signal_i]]
    b0=result[[iter]][[signal_i+4]]
    b_validation_test_lambda_R_nonlinear = result[[iter]][[signal_i+12]]
    y_hat= b0 + crossprod(matrix(tildePhiX_data, c(prod(dim(BB)), n_train)), as.vector(BB))
    res_mat[iter,signal_i] = sum((y_hat-y_all[[signal_i]][id_train])^2)
  }
  
  
  beta = svd_to_cp(BBprod[[2]][1,,])
  tildealpha0scale =  sqrt(sum(coef_truspline[[2]][-1]^2))
  tildealpha0 = coef_truspline[[2]][-1]/tildealpha0scale
  beta[[1]] = beta[[1]] * tildealpha0scale
  beta[[3]] = cbind(tildealpha0,tildealpha0,tildealpha0,tildealpha0)
  beta = rescale(beta, rescale_all = 0, rescale_l2 = 1, elastic = 1,alpha = 0.5)
  
  index_best = which.min(b_validation_test_lambda_R_nonlinear[,2])
  lambda_1 = b_validation_test_lambda_R_nonlinear[index_best, 6]
  lambda_2 = b_validation_test_lambda_R_nonlinear[index_best, 5]
  print(c(lambda_1,lambda_2))
  for(signal_i in 2){
    res_mat_app[iter,signal_i] = sum((y_all[[signal_i]][id_train] - yysplinepre[[signal_i]][id_train])^2)  +  lambda_1 * penaltyvalue(beta = beta, alpha = lambda_2, manifold = 0, QCQP = 1)
  }
  
}

res = list(res_mat=res_mat, res_mat_app=res_mat_app)
name = str_c("res",n_use, "_smaller_oj.Rdata")
save(res, file=name)
