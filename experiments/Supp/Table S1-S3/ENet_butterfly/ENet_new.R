#simulaion, enastic net
source('data_generation_different.R')
source('~/Desktop/Simulation/simulation_test/Rdata_to_MATLAB .R')




library(glmnet)




















y_all=matrix(0,4,1000)
y_all[,1:700]=y_train_MATLAB
y_all[,701:1000]=y_test_MATLAB
rm(y_test_MATLAB)
rm(y_train_MATLAB)

X_all=matrix(0,1000,4096)
X_all[1:700,1:4096]=X_data_train_MATLAB
X_all[701:1000,1:4096]=X_data_test_MATLAB
#rm(X_data_train_MATLAB)
#rm(X_data_test_MATLAB)
gc()


MSE_pre=matrix(0,4,50)
MSE_vali=c()

n_train=600
n_vali=150
lambda_seq=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)/n_train
alpha_seq=c(0,0.5,1)
#alpha_seq=c(0)
#alpha_seq=c(1)


num_lambda=length(lambda_seq)
num_alpha=length(alpha_seq)

num_tuning=num_lambda*num_alpha

index_use=c()

set.seed(2019)
id_matrix=matrix(0,100,1000)
for(iter in 1:100){
  id_matrix[iter,]=sample(1000,1000)  
}

result=list()
BB_idx=list()

for(iter in 1:50){
for(signal_i in 1:4){

  
  idtrain=id_matrix[iter,1:n_train]
  idtest=id_matrix[iter,(1000-n_vali+1):1000]
  
  X_train_dir=X_all[idtrain,1:4096]
  X_test_dir=X_all[idtest,1:4096]
  
  y_train_dir=y_all[signal_i,idtrain]
  y_test_dir=y_all[signal_i,idtest]
  

  
  library(glmnet)
  index_lambda=matrix(0,2,num_alpha)
  for(j in 1:num_alpha){
    for(i in 1:num_lambda){
      fit_dir<-glmnet(X_train_dir[1:n_train,],y_train_dir[1:n_train],family="gaussian",alpha=alpha_seq[j],lambda=lambda_seq[i],standardize = F)
      y_hat_dir_vali=X_test_dir[1:n_vali,]%*%as.matrix(coef(fit_dir)[-1])+coef(fit_dir)[1]
      MSE_vali[i]=sum((y_test_dir[(1:n_vali)]-y_hat_dir_vali)^2)/n_vali
      index1=which.min(MSE_vali)
      index_lambda[1,j]=index1
      index_lambda[2,j]=min(MSE_vali)
    }
    
  }
  index_2=which.min(index_lambda[2,])
  index_1=index_lambda[1,index_2]
  
  fit_dir<-glmnet(X_train_dir[1:n_train,],y_train_dir[1:n_train],family="gaussian",alpha=alpha_seq[index_2],lambda=lambda_seq[index_1],standardize = F)
  y_hat_dir_test=X_test_dir[1:n_vali,]%*%as.matrix(coef(fit_dir)[-1])+coef(fit_dir)[1]
  MSE_pre[signal_i,iter]=sum((y_test_dir[1:n_vali]-y_hat_dir_test)^2)/n_vali
  BB_idx[[signal_i]]=as.matrix(coef(fit_dir)[-1])
  BB_idx[[signal_i+4]]=coef(fit_dir)[1]
  BB_idx[[signal_i+8]]=MSE_pre[signal_i,iter]
  
}
  result[[iter]]= BB_idx
}

#save(MSE_pre,file='MSE_pre_ENet_700_a03.Rdata')
#save(result,file='ENet_500_a03.Rdata')
#save(MSE_pre,file='MSE_pre_ENet_400_a03.Rdata')
#id=4
#mean(MSE_pre[id,])
#sqrt(var(MSE_pre[id,]))

