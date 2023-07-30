validation_ENet <- function(alpha, lambda, X_train, y_train, X_vali, y_vali, X_test, y_test){
  d=dim(X_train)

  n=length(y_train)
  n_vali=length(y_vali)
  n_test=length(y_test)
  
  if(length(d)==3){
    
    X_train=sapply(c(1:n), function(i, XX) as.vector(XX[,,i]), XX=X_train)
    X_train=t(X_train)
    
    X_vali=sapply(c(1:n_vali), function(i, XX) as.vector(XX[,,i]), XX=X_vali)
    X_vali=t(X_vali)
    
    X_test=sapply(c(1:n_test), function(i, XX) as.vector(XX[,,i]), XX=X_test)
    X_test=t(X_test)
  
  }else if(length(d)==4){
    X_train=sapply(c(1:n), function(i, XX) as.vector(XX[,,,i]), XX=X_train)
    X_train=t(X_train)
    
    X_vali=sapply(c(1:n_vali), function(i, XX) as.vector(XX[,,,i]), XX=X_vali)
    X_vali=t(X_vali)
    
    X_test=sapply(c(1:n_test), function(i, XX) as.vector(XX[,,,i]), XX=X_test)
    X_test=t(X_test)
    
  }
  
 
  lambda=lambda/n
  
  num_lambda=length(lambda)
  num_alpha=length(alpha)
  num_all=num_lambda*num_alpha
  

  
  
  
  b_validation_test_lambda_R_nonlinear=matrix(0,num_all,5)
  BBbig_nonlinear=list()
  

    for(i2 in 1:num_alpha){
      for(i3 in 1:num_lambda){  
        fit <- glmnet(X_train, y_train ,family="gaussian",alpha=alpha[i2],lambda=lambda[i3],standardize = F)
        yhat_vali <- X_vali %*% as.matrix(coef(fit)[-1]) + coef(fit)[1]
        yhat_test <- X_test %*% as.matrix(coef(fit)[-1]) + coef(fit)[1]
        MSE_vali <- sum((yhat_vali - y_vali)^2)/n_vali
        MSE_test <- sum((yhat_test - y_test)^2)/n_test
        
        num_location=(i2-1)*num_lambda+i3
        
        BBbig_nonlinear[[num_location]]=as.matrix(coef(fit)[-1])
        b_validation_test_lambda_R_nonlinear[num_location,1]=coef(fit)[1]
        b_validation_test_lambda_R_nonlinear[num_location,2]=MSE_vali
        b_validation_test_lambda_R_nonlinear[num_location,3]=MSE_test
       
        b_validation_test_lambda_R_nonlinear[num_location,4]=alpha[i2]
        b_validation_test_lambda_R_nonlinear[num_location,5]=lambda[i3]
        
        
      }
    }
    
  index_best=which.min(b_validation_test_lambda_R_nonlinear[,2])
  index_best=index_best[1]
  
  b_final=b_validation_test_lambda_R_nonlinear[index_best,1]
  BB_final=BBbig_nonlinear[[index_best]]
  MSE_vali_final=b_validation_test_lambda_R_nonlinear[index_best,2]
  MSE_pre_final=b_validation_test_lambda_R_nonlinear[index_best,3]


  
  res <- list(bb=b_final, BB=BB_final, MSE_vali=MSE_vali_final, MSE_pre=MSE_pre_final, b_validation_test_lambda_R_nonlinear=b_validation_test_lambda_R_nonlinear)
  
  return(res)
  
  
}

