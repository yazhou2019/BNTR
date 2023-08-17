library(stringr)

ndata = 1000
distri_seq = c("uniform", 'tnorm')
distri = distri_seq[1]
isbutterfly = 1



# load the data set
if(isbutterfly == 0){
  loadname = str_c("./SimResults/ENet_", ndata, str_c("_linear_new_", distri,".Rdata"))
}else if(isbutterfly == 1){
  # 默认就是原始的uniform setting的结果
  loadname = str_c("./SimResults/ENet_", ndata, str_c("_linear_butterfly",".Rdata"))
}

load(loadname)


source('./ComponentsLin/functions_needed.R')
source('./ComponentsLin/validation_result_linear.R')
source('./ComponentsLin/broadcasted_sparsetenreg_linear.R')




if(isbutterfly == 0){
  if(distri=="uniform"){
    source('./SimDataGeneration/true_signal_different.R')
    source('./SimDataGeneration/data_generation_different.R')
  }
  if(distri=="tnorm"){
    source('./SimDataGeneration/true_signal_different_tnorm.R')
    source('./SimDataGeneration/data_generation_different_tnorm.R')  
  }
}else if(isbutterfly == 1){
  source('./SimDataGeneration/true_signal_different_butterfly.R')
  source('./SimDataGeneration/data_generation_different_butterfly.R')
}

rm(X_data_SNR)
rm(X_datacross_reg_SNR)
gc()


set.seed(20210000)
n_MC = 10000



if(distri=="uniform"){
  X_data=array(runif(64*64*n_MC,0,1),c(64,64,n_MC))
}

if(distri=="tnorm"){
  X_data = array(rtnorm(64*64*n_MC, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,n_MC))
  #knots = quantile(c(X_data), probs = c(seq(0, 1, 1/(num_knots - 1))))
  #X_data=array(runif(64*64*n_MC,0,1),c(64,64,n_MC))
}


if(isbutterfly == 0){
  X_datacross_reg_test=list()
  X_datacross_reg_test[[1]]=X_data
  X_datacross_reg_test[[2]]=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
  X_datacross_reg_test[[3]]=X_data+adjust_nonlinear*0.5*cos(2*pi*X_data)
  
  y_true=list()
  y_true[[1]]=1+t(ctprod(BBprod[[1]],X_datacross_reg_test[[1]],2))
  y_true[[2]]=1+t(ctprod(BBprod[[2]],X_datacross_reg_test[[2]],2))
  y_true[[3]]=1+t(ctprod(BBprod[[3]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[4]],X_datacross_reg_test[[2]],2))
  y_true[[4]]=1+t(ctprod(BBprod[[5]],X_datacross_reg_test[[2]],2))+t(ctprod(BBprod[[6]],X_datacross_reg_test[[3]],2))
  
  if(distri=="tnorm"){
    BBprod_butterfly =array(0,c(1,dim(butterfly)))
    BBprod_butterfly[1,,] =  butterfly
    y_true[[5]]=1+t(ctprod(BBprod_butterfly,X_datacross_reg_test[[2]],2))
  }
}else if(isbutterfly ==1){
  X_datacross_reg_test=X_data+adjust_nonlinear*sin(2*pi*(X_data-0.5)^2)
  y_true=1+t(ctprod(BBprod,X_datacross_reg_test,2))
}


#y_hat = list()


res_mat = matrix(0, 50, 5)
for(iter in 1:50){
  
  if(isbutterfly==0){  
    
    if(distri == "uniform"){
      for(signal_i in 1:4){
        BB=result[[iter]][[signal_i]]
        b0=result[[iter]][[signal_i+4]]
        BBprod = array(0, c(1,64,64))
        BBprod[1,,] = BB
        y_hat = b0+t(ctprod(BBprod,X_data,2))
        res_mat[iter,signal_i] = sum((y_hat-y_true[[signal_i]])^2)/n_MC
      }
    }
    
    if(distri == "tnorm"){  
      for(signal_i in 1:5){
        BB=result[[iter]][[signal_i]]
        b0=result[[iter]][[signal_i+5]]
        BBprod = array(0, c(1,64,64))
        BBprod[1,,] = BB
        y_hat = b0+t(ctprod(BBprod,X_data,2))
        res_mat[iter,signal_i] = sum((y_hat-y_true[[signal_i]])^2)/n_MC
      }
    }
  }else if(isbutterfly==1){
    BB=result[[iter]][[1]]
    b0=result[[iter]][[1+4]]
    #y_hat= b0 + crossprod(matrix(tildePhiX_data, c(prod(dim(BB)), n_MC)), as.vector(BB))
    BBprod = array(0, c(1,64,64))
    BBprod[1,,] = BB
    y_hat = b0+t(ctprod(BBprod,X_data,2))
    res_mat[iter,5] = sum((y_hat-y_true)^2)/n_MC
  }
  
  
}



res_final = matrix(0, 2, 5)
for(signal_i in 1:5){
  res_final[1,signal_i] = mean(res_mat[,signal_i])
  res_final[2,signal_i] = sqrt(var(res_mat[,signal_i]))
}  
res_df <- data.frame(round(res_final,4))
colnames(res_df) <- c('case 1', 'case 2','case 3','case 4','butterfly')
rownames(res_df) <- c('mean', 'std')
print(res_df)

if(isbutterfly==0){
  savename= str_c('ENet_res_mat', ndata, str_c("_linear_new_",distri,".Rdata"))
}else if(isbutterfly==1){
  savename= str_c('ENet_res_mat', ndata, str_c("_linear_butterfly_",".Rdata"))
}


save(file=savename, res_mat)
