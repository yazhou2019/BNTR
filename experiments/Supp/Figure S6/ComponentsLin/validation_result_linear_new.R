library(stringr)
validation_result<-function(R,alpha,lambda,X_data_train,y_train,X_data_vali,y_vali,X_data_test,y_test){

  n=length(y_train)
  lambda=lambda/n

  
num_R=length(R)
num_lambda=length(lambda)
num_alpha=length(alpha)
num_all=num_R*num_lambda*num_alpha

b_validation_test_lambda_R_nonlinear=matrix(0,num_all,6)
BBbig_nonlinear=list()
  

  
for(i1 in 1:num_R){
for(i2 in 1:num_alpha){
for(i3 in 1:num_lambda){  

if(i3==1){
  
  fit=try(quiet(broadcasted_sparsetenreg(TolFun=0.0001,rescale_all=1,rescale_l2=1,r=R[i1],X_data_train,y_train,penalty=c("L1"),knots=NA,restriction = 0,lambda=lambda[i3],gamma=lambda[i3],alpha_gamma=alpha[i2],
                                         alpha=alpha[i2],Replicates = 5,family="gaussian",MaxIter = 10000,manifold=0,warmstart=0, QCQP=0)))
  if("try-error" %in% class(fit)){
  }else{
    fit_nonlinear=fit
  }
  
}else{
  fit=try(quiet(broadcasted_sparsetenreg(TolFun=0.0001,rescale_all=1,rescale_l2=1,r=R[i1],X_data_train,y_train,lambda=lambda[i3],gamma=lambda[i3],alpha_gamma=alpha[i2],
                                                                              alpha=alpha[i2],Replicates = 5,family="gaussian",MaxIter = 10000,manifold=0,warmstart=0, QCQP=0,B0=fit_nonlinear$beta,beta0 =fit_nonlinear$beta0)))
  
  if("try-error" %in% class(fit)){
  }else{
    fit_nonlinear=fit
  }
}

num_location=(i1-1)*num_alpha*num_lambda+(i2-1)*num_lambda+i3
printout3=str_c('Complete ','%',round(100*num_location/num_all,2))
print(printout3)
  
  
  
b0=fit_nonlinear$beta0
BB=full_R(fit_nonlinear$beta)
#MSE_nonlinear_validation=prediction_function_linear(b0,BB,X_data_vali,y_vali,family='gaussian')
#MSE_nonlinear_test=prediction_function_linear(b0,BB,X_data_test,y_test,family='gaussian')


MSE_nonlinear_validation=prediction_function_linear_new(b0,BB,X_data_vali,y_vali,family='gaussian')
MSE_nonlinear_test=prediction_function_linear_new(b0,BB,X_data_test,y_test,family='gaussian')


BBbig_nonlinear[[num_location]]=BB
b_validation_test_lambda_R_nonlinear[num_location,1]=b0
b_validation_test_lambda_R_nonlinear[num_location,2]=MSE_nonlinear_validation
b_validation_test_lambda_R_nonlinear[num_location,3]=MSE_nonlinear_test
b_validation_test_lambda_R_nonlinear[num_location,4]=R[i1] 
b_validation_test_lambda_R_nonlinear[num_location,5]=alpha[i2]
b_validation_test_lambda_R_nonlinear[num_location,6]=lambda[i3]

}  
}
}

#avoid numerical problems
for(num_location in 1:num_all){
  if(is.nan(b_validation_test_lambda_R_nonlinear[num_location,2])==1){
    b_validation_test_lambda_R_nonlinear[num_location,2]=Inf  
  }
  if(is.na(b_validation_test_lambda_R_nonlinear[num_location,2])==1){
    b_validation_test_lambda_R_nonlinear[num_location,2]=Inf  
  }
  if(b_validation_test_lambda_R_nonlinear[num_location,2]==0){
    b_validation_test_lambda_R_nonlinear[num_location,2]=Inf   
  }
  
}

index_best=which.max(b_validation_test_lambda_R_nonlinear[,2])

res=list(BB=BBbig_nonlinear[[index_best]],MSE_vali=b_validation_test_lambda_R_nonlinear[index_best,2],MSE_pre=b_validation_test_lambda_R_nonlinear[index_best,3],knots=knots,BBbig_nonlinear=BBbig_nonlinear,b_validation_test_lambda_R_nonlinear=b_validation_test_lambda_R_nonlinear,bb=b_validation_test_lambda_R_nonlinear[index_best,1])  
  
return(res)
}

quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 
