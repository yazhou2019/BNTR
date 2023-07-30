library(glmnet)
source('./ENetComponents/ENet_source.R')



BB_idx=list()
result=list()

for(idx in 1:50){


#all the data set
set.seed(idx*100)


load("./SimResults/Simul_n1000_rep50_final_fix_new.Rdata")
BB_list = data_all[[51]]
v_list = data_all[[52]]
data_all = data_all[[idx]]
X_data = data_all$X
y_all = data_all$y


for(signal_i in 1:5){
  
  id_train=1:n_train
  id_vali=(n-n_vali+1):n
  id_test=851:1000
  
  res=validation_ENet(alpha,lambda,X_data[,,id_train],y_all[[signal_i]][id_train],X_data[,,id_vali],y_all[[signal_i]][id_vali],X_data[,,id_test],y_all[[signal_i]][id_test])
  #BB_idx[[signal_i]]=matrix(res$BB, 64,64)
  #BB_idx[[signal_i+4]]=res$bb
  #BB_idx[[signal_i+8]]=res$MSE_pre
  
  BB_idx[[signal_i]]=matrix(res$BB, 64,64)
  BB_idx[[signal_i+5]]=res$bb
  BB_idx[[signal_i+10]]=res$MSE_pre
  BB_idx[[signal_i+15]] = "ENet"
  BB_idx[[signal_i+20]] = BB_list[[signal_i]]
  BB_idx[[signal_i+25]] = v_list[[signal_i]]
  
  

}


result[[idx]]=BB_idx

}

filename <- str_c("ENet_", n_use, "_new_20230302.Rdata")
setwd("./SimResults")
save(result,file = filename)
setwd("../")


