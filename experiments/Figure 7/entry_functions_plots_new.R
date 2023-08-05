



# TLR-rescaled
# BB_idx[[1]] = res$MSE_vali
# BB_idx[[2]] = res$MSE_pre
# BB_idx[[3]] = res$BB
# BB_idx[[4]] = res$bb
# BB_idx[[5]] = res$knots


# ENEt
# BB_idx[[1]]=matrix(res$BB, 64,10,10)
# BB_idx[[1+5]]=res$bb
# BB_idx[[1+10]]=res$MSE_pre
# BB_idx[[1+15]] = "ENet"

# BNTR
# BB_idx[[1]] = res$MSE_vali
# BB_idx[[2]]=res$MSE_pre
# BB_idx[[3]]=res$b_validation_test_lambda_R_nonlinear
# index_best = which.min(res$b_validation_test_lambda_R_nonlinear[,2])
# BB_idx[[4]]=res$betabig[[index_best]]

#source("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/ComponentsNonLin/functions_needed.R")

source("./ComponentsNonLin/functions_needed.R")

obtain_median_BNTR_iter <-function(){
  #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/data/monkey_1_new_K11_tuning.Rdata")
  load("./data/monkey_1_new_K11_tuning.Rdata")
  mse_mat = matrix(0,10,2)
  for(iter in 1:10){
    mse_mat[iter, 1] = result[[iter]][[1]]
    mse_mat[iter, 2] = result[[iter]][[2]]
  }
  iter_index = which(order(mse_mat[,2])==5)
  return(iter_index)
}

sum_axis_keepdim <-function(BB, tildePhiX_i){
  num_basis = dim(BB)[length(dim(BB))]
  res = array(0, c(dim(BB))[-length(dim(BB))])
  for(i in 1:num_basis){
    res = res + tildePhiX_i[,,,i] * BB[,,,i]
  }
  return(res)
}



obtain_frue <- function(){
  X_seq = seq(-2.75,2.85,0.05)
  X = array(0, c(64,10,10, length(X_seq)))
  y = array(0, c(64,10,10, length(X_seq)))
  for(i1 in 1:64){
    for(i2 in 1:10){
      for(i3 in 1:10){
        X[i1,i2,i3, ] = X_seq
      }
    }
  }
  #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/data/truepara.Rdata")
  load("./data/truepara.Rdata")
  #bb = truepara$b
  BB = truepara$BB
  #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/data/knots_used.Rdata")
  load("./data/knots_used.Rdata")
  tildePhiX = tildePhiX_trans(X, knots_used, order=4)
  for(i in 1:length(X_seq)){
    y[,,,i] = sum_axis_keepdim(BB, tildePhiX[,,,,i]) #+ bb
  }
  
  y_mean = array(0, c(64,10,10))
  for(i1 in 1:64){
    for(i2 in 1:10){
      for(i3 in 1:10){
        y_mean[i1,i2,i3] =  mean(y[i1,i2,i3,])
        y[i1,i2,i3,] =  y[i1,i2,i3,] - y_mean[i1,i2,i3]
      }
    }
  }
  
  
  return(y)
}


obtain_fhat <- function(method="BNTR", iter=1){
  X_seq = seq(-2.75,2.85,0.05)
  X = array(0, c(64,10,10, length(X_seq)))
  y = array(0, c(64,10,10, length(X_seq)))
  for(i1 in 1:64){
    for(i2 in 1:10){
      for(i3 in 1:10){
        X[i1,i2,i3, ] = X_seq
      }
    }
  }
  
  if(method=="BNTR"){
    #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/data/knots_used_list.Rdata")
    #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/monkey_1_new_K11_tuning_sim_new2.Rdata")
    load("./data/knots_used_list.Rdata")
    load("./SimResults/monkey_1_new_K11_tuning_sim_new2.Rdata")
    bb = result[[iter]][[3]][which.min(result[[iter]][[3]][,2]),1]
    BB = full_R(result[[iter]][[4]])
    tildePhiX = tildePhiX_trans(X, knots_used[[iter]], order=4)
    for(i in 1:length(X_seq)){
      y[,,,i] = sum_axis_keepdim(BB, tildePhiX[,,,,i]) + bb
    }
  }
  
  if(method=="TLR_rescaled"){
    #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/TLR2_linear_new_sim_new2.Rdata")
    load("./SimResults/TLR2_linear_new_sim_new2.Rdata")
    BB = result[[iter]][[3]]
    for(i in 1:length(X_seq)){
      y[,,,i]= BB * X[,,,i]
    }
    
  }
  
  if(method=="ENet"){
    #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/ENet_real_20230103_2.Rdata")
    load("./SimResults/ENet_real_20230103_2.Rdata")
    BB = result[[iter]][[iter]][[1]]
    for(i in 1:length(X_seq)){
      y[,,,i]= BB * X[,,,i]
    }
  }
  
  if(method=="TLR"){
    library(R.matlab)
    #tem = readMat("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/BB_all_RealDataSim2.mat")
    tem = readMat("./SimResults/BB_all_RealDataSim2.mat")
    BB=  tem$BB.all[,,,iter]
    for(i in 1:length(X_seq)){
      y[,,,i]= BB * X[,,,i]
    }
    
  }
  
  
  
  y_mean = array(0, c(64,10,10))
  for(i1 in 1:64){
    for(i2 in 1:10){
      for(i3 in 1:10){
        y_mean[i1,i2,i3] =  mean(y[i1,i2,i3,])
        y[i1,i2,i3,] =  y[i1,i2,i3,] - y_mean[i1,i2,i3]
      }
    }
  }
  return(y)
}


obtain_entry_l2norm <-function(yhat_ten_list, ytrue_ten){
  entry_l2norm = array(0, c(64,10,10,10))
  for(iter in 1:10){
    for(i1 in 1:64){
      for(i2 in 1:10){
        for(i3 in 1:10){
          entry_l2norm[i1,i2,i3, iter] =sum(abs(yhat_ten_list[[iter]][i1,i2,i3, ] - ytrue_ten[i1,i2,i3,]))
        }
      }
    }
    
  }
  return(entry_l2norm)
}



obtain_BB <-function(method, iter){
  
  if(method=="TLR_rescaled"){
    #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/TLR2_linear_new_sim_new2.Rdata")
    load("./SimResults/TLR2_linear_new_sim_new2.Rdata")
    BB = result[[iter]][[3]]
  }
  
  if(method=="ENet"){
    #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/ENet_real_20230103_2.Rdata")
    load("./SimResults/ENet_real_20230103_2.Rdata")
    BB = result[[iter]][[iter]][[1]]
  }
  
  if(method=="TLR"){
    library(R.matlab)
    #tem = readMat("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/BB_all_RealDataSim2.mat")
    tem = readMat("./SimResults/BB_all_RealDataSim2.mat")
    BB=  tem$BB.all[,,,iter]
  }
  
  
  return(BB)
  
}



obtain_index_list <- function(entry_l2norm, iter_index, qnum=c(0.5, 0.75, 1), nonzero=FALSE, method_name=NA){
  
  if(nonzero==TRUE){
    BB = obtain_BB(method=method_name, iter=iter_index)
    BB_vec = c(BB) 
    BB_vec_nonzero = BB_vec[BB_vec!=0]
    nonzeroindex = which(BB_vec!=0)
    
    entry_l2norm_vec = c(entry_l2norm[,,,iter_index])
    entry_l2norm_vec_nonzero = entry_l2norm_vec[nonzeroindex]
    
    temp0 = round(quantile(1:length(nonzeroindex), qnum))
    temp = c()
    for(i in 1:length(temp0)){
      temp[i] = nonzeroindex[order(entry_l2norm_vec_nonzero)[temp0[i]]]
    }
  }else{
    temp0 = round(quantile(1:6400, qnum))
    temp = c()
    for(i in 1:length(temp0)){
      temp[i] = order(entry_l2norm[,,,iter_index])[temp0[i]]
    }
  }
  
  index_list = list()
  temp2 = array(1:6400, c(64,10,10))
  for(i in 1:length(temp)){
    for(i1 in 1:64){
      for(i2 in 1:10){
        for(i3 in 1:10){
          if(temp2[i1,i2,i3]==temp[i]){
            index_list[[i]] = c(i1,i2,i3)
            print(c(entry_l2norm[i1,i2,i3,iter_index], temp0[i]))
          }
        }
      }
    }
  }
  
  return(index_list)
  
}


obtain_index_list_avg <- function(entry_l2norm, qnum=c(0.5, 0.75, 1)){
  entry_l2norm2 = array(0, c(64,10,10))
  for(i1 in 1:64){
    for(i2 in 1:10){
      for(i3 in 1:10){
        entry_l2norm2[i1,i2,i3] = mean(entry_l2norm[i1,i2,i3,])
      }
    }
  }
  
  
  
  temp0 = round(quantile(1:6400, qnum))
  temp = c()
  for(i in 1:length(temp0)){
    temp[i] = order(c(entry_l2norm2))[temp0[i]]  
  }
  
  index_list = list()
  temp2 = array(1:6400, c(64,10,10))
  for(i in 1:length(temp)){
    for(i1 in 1:64){
      for(i2 in 1:10){
        for(i3 in 1:10){
          if(temp2[i1,i2,i3]==temp[i]){
            index_list[[i]] = c(i1,i2,i3)
            print(c(entry_l2norm2[i1,i2,i3],temp0[i]))
          }
        }
      }
    }
  }
  
  return(index_list)
  
}

obtain_index_list_mean_iter <- function(entry_l2norm, qnum=0.55){
  index_list = list()
  for(iter_index in 1:10){
    temp = round(quantile(1:6400, qnum))
    temp = order(entry_l2norm[,,,iter_index])[temp]  
    temp2 = array(1:6400, c(64,10,10))
    for(i1 in 1:64){
      for(i2 in 1:10){
        for(i3 in 1:10){
          if(temp2[i1,i2,i3]==temp){
            index_list[[iter_index]] = c(i1,i2,i3)
          }
        }
      }
    }
  }
  return(index_list)
}

obtain_index_list_BB <- function(qnum=1, abs_o=FALSE, method_name="TLR_rescaled"){
  index_list = list()
  for(iter_index in 1:10){
    temp = round(quantile(1:6400, qnum))
    if(method_name == "TLR_rescaled"){
      #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/TLR2_linear_new_sim_new2.Rdata")
      load("./SimResults/TLR2_linear_new_sim_new2.Rdata")
      BB = result[[iter_index]][[3]]
    }
    if(method_name == "ENet"){
      #load("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim/SimResults/ENet_real_20230103_2.Rdata")
      load("./SimResults/ENet_real_20230103_2.Rdata")
      BB = result[[iter_index]][[iter_index]][[1]]  
    }
    
    if(abs_o ==TRUE){
      BB = abs(BB)
    }
    temp = order(BB)[temp]  
    temp2 = array(1:6400, c(64,10,10))
    for(i1 in 1:64){
      for(i2 in 1:10){
        for(i3 in 1:10){
          if(temp2[i1,i2,i3]==temp){
            index_list[[iter_index]] = c(i1,i2,i3)
          }
        }
      }
    }
  }
  return(index_list)
}

library(ggplot2)


obtain_ggplot_BNTR_new <- function(used_index, ytrue_ten, yhat_ten_list, method_name="BroadcasTR", lr=c(0,0)){
  #used_index = c(44,5,1)
  
  method_name = "Estimation"
  method_names = c(method_name)
  #size_names = c("n=1000")
  small_index_tem = 1:113
  small_index = small_index_tem
  
  res_signal = list()
  
  all_value = c()
  
  ob_value = ytrue_ten[used_index[1],used_index[2],used_index[3],]
  all_value = c(all_value, ob_value)
  iterations = seq(-2.75,2.85,0.05)
  num_xpoints = 113
  groups = c(rep('True', num_xpoints))
  df_matrix = matrix(0, length(small_index),15) 
  df_matrix[,1] = iterations[small_index]
  df_matrix[,2] = ob_value[small_index]
  
  df_list = list()
  df_list[[1]]=iterations[small_index]
  df_list[[2]] = ob_value[small_index]
  
  for(iter_i in 1:10){
    df_list[[iter_i + 2]] = yhat_ten_list[[iter_i]][used_index[1],used_index[2],used_index[3],small_index]
    all_value = c(all_value, df_list[[iter_i + 2]])
    df_matrix[, iter_i + 2] = yhat_ten_list[[iter_i]][used_index[1],used_index[2],used_index[3],small_index]
  }  
  
  df_matrix[,14] = 1
  df_matrix[,15] = 2
  
  df_small = data.frame(df_matrix)
  df_small$X14=rep("Simulated Truth",1)
  df_small$X15=rep(method_name,1)
  
  
  
  df = df_small
  #res = ggplot(df,aes(x=iteration, y=ob_value,colour=group, linetype=group, group=group))+geom_line(lwd = 0.5) + labs(x=NULL, y="") #+labs(x=NULL,y=str_c("Case ", signal_i)) 
  alpha_value = 0.2
  #res = ggplot(df,aes(x=eval(parse(text='X1')), y=eval(parse(text='X2')),group=eval(parse(text='X54')),linetype=eval(parse(text='X54')),color=eval(parse(text='X54'))))+geom_line(lwd = 0.5) + labs(x=NULL, y="") 
  res = ggplot(df,aes(x=X1, y=X2,group=factor(X14),linetype=X14,color=X14))+ labs(x=NULL, y="") #+geom_line(lwd = 0.5) 
  
  res = res + geom_line(data=df,aes(y=X3,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X4,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X5,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X6,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X7,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X8,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X9,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X10,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X11,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X12,group=X15,linetype=X15,color=X15),alpha=alpha_value)
  res = res + geom_line(data=df,aes(y=X2,group=X14,linetype=X14,color=X14),alpha=1,lwd = 0.5)
  res = res + scale_color_manual(name="Method",values=c("red","black"))
  
  res = res + theme(legend.title = element_text(size=16),legend.text = element_text(size=12))
  res = res + scale_linetype_manual(name="Method",values = c("dashed","solid")) 
  if (sum(lr)==0){
    res = res + ylim(min(all_value)-0.05,max(all_value)+0.05) #+geom_point() #+ labs(title = size_names[size_i]) + theme(plot.title=element_text(hjust=0.5))
  }else{
    res = res + ylim(lr[0],lr[1])  
  }
  return(res)
  
}




index_list_func <- function(method_name, qnum,nonzero=FALSE){
  ytrue_ten = obtain_frue()
  yhat_ten_list = list()
  for(i in 1:10){
    yhat_ten_list[[i]] = obtain_fhat(method=method_name, iter=i)
  }
  entry_l2norm = obtain_entry_l2norm(yhat_ten_list, ytrue_ten)
  iter_index = obtain_median_BNTR_iter()
  
  index_list = obtain_index_list(entry_l2norm, iter_index, qnum = c(qnum), nonzero = nonzero, method_name = method_name)
  return(index_list)
}

index_list_func_avg <- function(method_name, qnum){
  ytrue_ten = obtain_frue()
  yhat_ten_list = list()
  for(i in 1:10){
    yhat_ten_list[[i]] = obtain_fhat(method=method_name, iter=i)
  }
  entry_l2norm = obtain_entry_l2norm(yhat_ten_list, ytrue_ten)
  
  iter_index = obtain_median_BNTR_iter()
  
  index_list = obtain_index_list_avg(entry_l2norm,qnum = qnum)
  return(index_list)
}



obtain_ggplot <-function(index_list, method_name, lr=c(0,0)){
  ytrue_ten = obtain_frue()
  yhat_ten_list = list()
  for(i in 1:10){
    yhat_ten_list[[i]] = obtain_fhat(method=method_name, iter=i)
  }
  res = list()
  for(i in 1:length(index_list)){
    res[[i]] = obtain_ggplot_BNTR_new(index_list[[i]], ytrue_ten, yhat_ten_list, method_name = method_name, lr=lr)
  }
  return(res)
}



#index_list = index_list_func("BNTR", 0.5555) 
#index_list = index_list_func("ENet", 0.5555) 

#index_list = index_list_func_avg("TLR_rescaled", 0.6) 
#index_list = index_list_func_avg("BNTR", 0.6) 

# 0.6 work, 0.5555 work
res = list()
index_list = index_list_func("TLR", 0.6, nonzero = TRUE) 
res[[1]] = obtain_ggplot(index_list, "TLR")[[1]] + ylim(-0.1, 0.3)
res[[2]] = obtain_ggplot(index_list, "TLR_rescaled")[[1]] + ylim(-0.1, 0.3)
res[[3]] = obtain_ggplot(index_list, "ENet")[[1]] + ylim(-0.1, 0.3)
res[[4]] = obtain_ggplot(index_list, "BNTR")[[1]] + ylim(-0.1, 0.3)

index_list = index_list_func("TLR_rescaled", 0.6,nonzero = TRUE) 
res[[5]] = obtain_ggplot(index_list, "TLR")[[1]] + ylim(-0.1, 0.3)
res[[6]] = obtain_ggplot(index_list, "TLR_rescaled")[[1]] + ylim(-0.1, 0.3)
res[[7]] = obtain_ggplot(index_list, "ENet")[[1]] + ylim(-0.1, 0.3)
res[[8]] = obtain_ggplot(index_list, "BNTR")[[1]] + ylim(-0.1, 0.3)


index_list = index_list_func("ENet", 0.6,nonzero = TRUE) 
res[[9]] = obtain_ggplot(index_list, "TLR")[[1]] + ylim(-0.1, 0.3)
res[[10]] = obtain_ggplot(index_list, "TLR_rescaled")[[1]] + ylim(-0.1, 0.3)
res[[11]] = obtain_ggplot(index_list, "ENet")[[1]] + ylim(-0.1, 0.3)
res[[12]] = obtain_ggplot(index_list, "BNTR")[[1]] + ylim(-0.1, 0.3)


index_list = index_list_func("BNTR", 0.6) 
res[[13]] = obtain_ggplot(index_list, "TLR")[[1]] + ylim(-0.1, 0.3)
res[[14]] = obtain_ggplot(index_list, "TLR_rescaled")[[1]] + ylim(-0.1, 0.3)
res[[15]] = obtain_ggplot(index_list, "ENet")[[1]] + ylim(-0.1, 0.3)
res[[16]] = obtain_ggplot(index_list, "BNTR")[[1]] + ylim(-0.1, 0.3)

# library(ggpubr)
# library(stringi)
# library(stringr)
# ggarrange(res[[1]] + labs(title = "TLR") + theme(plot.title=element_text(hjust=0.5,margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
#           res[[2]] + labs(title = "TLR-rescale")  + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
#           res[[3]] + labs(title = "ENet")+ theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
#           res[[4]] + labs(title = "BroadcasTR")  + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
#           res[[5]],res[[6]],res[[7]],res[[8]], 
#           res[[9]],res[[10]],res[[11]],res[[12]], 
#           res[[13]],res[[14]],res[[15]],res[[16]], 
#           ncol=4, nrow=4, common.legend = TRUE, legend='right')

#index_list = index_list_func("BNTR", 0.55) 
library(ggpubr)
library(stringi)
library(stringr)
#png("temtem.png",units="in", width=9.52/1.3, height=11/2,res=300)
ggarrange(res[[13]] + labs(title = "TLR") + theme(plot.title=element_text(hjust=0.5,margin = margin(t = 0, r = 0, b = 2.5, l = 0, unit = "pt"))),
          res[[14]] + labs(title = "TLR-rescale")  + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0 , r = 0, b = 2.5, l = 0, unit = "pt"))),
          res[[15]] + labs(title = "ENet")+ theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 2.5, l = 0, unit = "pt"))),
          res[[16]] + labs(title = "BroadcasTR")  + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 2.5, l = 0, unit = "pt"))),
          ncol=2, nrow=2, common.legend = TRUE, legend='right')
#dev.off()
















