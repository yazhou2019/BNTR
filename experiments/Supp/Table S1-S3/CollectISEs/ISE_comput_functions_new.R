#integration of the 
library(rTensor)

integration_of_entry_fun=function(BB,knots){
  
  p=dim(BB)
  d=length(p)
  BB_mean=array(0,c(p[1:(d-1)]))
  if(d==4){
    for(i in 1:p[1]){
      for(j  in 1:p[2]){
        for(k in 1:p[3]){
          coef=BB[i,j,k,]
          BB_mean[i,j,k]=coef[1]/2*(knots[5]^2-knots[1]^2)+coef[2]/3*(knots[5]^3-knots[1]^3)+coef[3]/4*(knots[5]^4-knots[1]^4)+coef[4]/4*(knots[5]-knots[1])^4+coef[5]/4*(knots[5]-knots[2])^4+coef[6]/4*(knots[5]-knots[3])^4
          
        }
      }
    }
  }
  
  if(d==3){
    for(i in 1:p[1]){
      for(j  in 1:p[2]){
        coef=BB[i,j,]
        BB_mean[i,j]=coef[1]/2*(knots[5]^2-knots[1]^2)+coef[2]/3*(knots[5]^3-knots[1]^3)+coef[3]/4*(knots[5]^4-knots[1]^4)+coef[4]/4*(knots[5]-knots[1])^4+coef[5]/4*(knots[5]-knots[2])^4+coef[6]/4*(knots[5]-knots[3])^4
      }
    }
  }
  
  return(BB_mean)
}


transform2D=function(BBB){
  p=dim(BBB)
  d=length(p)
  BB_l2=array(0,c(p[1:(d-1)]))
  
  if(d==3){
    for(i in 1:p[1]){
      for(j  in 1:p[2]){
        BB_l2[i,j]=sqrt(sum(BBB[i,j,]^2))
        
      }
    }
  }
  return(BB_l2)
  
}



sum_mean_true_signal=function(accuracy=0.00010001, adjust_nonlinear=0.6){
  sequse=seq(0,1,accuracy)
  X_data_SNR=array(0,c(64,64,length(sequse)))
  for(i in 1:64){
    for(j in 1:64){
      X_data_SNR[i,j,]=sequse
    }
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
  
  summean=list()
  
  summean[[1]]=mean(y_SNR[[1]]-1)
  summean[[2]]=mean(y_SNR[[2]]-1)
  summean[[3]]=mean(y_SNR[[3]]-1)
  summean[[4]]=mean(y_SNR[[4]]-1)
  
  localmean=list()
  localmean[[1]]=mean(X_datacross_reg_SNR[[1]])
  localmean[[2]]=mean(X_datacross_reg_SNR[[2]])
  localmean[[3]]=mean(X_datacross_reg_SNR[[3]])
  
  test=list(summean=summean,localmean=localmean)
  
  return(test)
}


true_entryfunction=function(accuracy=0.001,localmean){
  
  #accuracy=0.00010001
  
  
  source('./SimDataGeneration/true_signal_different.R')
  
  #noise level
  sigma=0.1
  adjust_nonlinear=0.3
  #presentation
  adjust_nonlinear=adjust_nonlinear*2
  
  
  sigma_use=c()
  
  sequse=seq(0,1,accuracy)
  X_data_SNR=array(0,c(64,64,length(sequse)))
  for(i in 1:64){
    for(j in 1:64){
      X_data_SNR[i,j,]=sequse
    }
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
  
  entry_function_matrix_list=list()
  entry_function_array=array(0,c(64,64,length(sequse)))
  
  for(signal_i_use in 1:4){
    print(signal_i_use)
    if(signal_i_use==1){
      for(i in 1:length(sequse)){  
        entry_function_array[,,i]=hadamard.prod(BBprod[[1]][1,,],X_datacross_reg_SNR[[1]][,,i])-hadamard.prod(BBprod[[1]][1,,],matrix(localmean[[1]],64,64))
      }
      entry_function_matrix_list[[signal_i_use]]=entry_function_array
    }
    if(signal_i_use==2){
      for(i in 1:length(sequse)){
        entry_function_array[,,i]=hadamard.prod(BBprod[[2]][1,,],X_datacross_reg_SNR[[2]][,,i])
        entry_function_array[,,i]=entry_function_array[,,i]-hadamard.prod(BBprod[[2]][1,,],matrix(localmean[[2]],64,64))
      }
      entry_function_matrix_list[[signal_i_use]]=entry_function_array
    }
    
    if(signal_i_use==3){
      for(i in 1:length(sequse)){
        entry_function_array[,,i]=hadamard.prod(BBprod[[3]][1,,],X_datacross_reg_SNR[[2]][,,i])+hadamard.prod(BBprod[[4]][1,,],X_datacross_reg_SNR[[2]][,,i])
        entry_function_array[,,i]=entry_function_array[,,i]-hadamard.prod(BBprod[[3]][1,,]+BBprod[[4]][1,,],matrix(localmean[[2]],64,64))
      }
      entry_function_matrix_list[[signal_i_use]]=entry_function_array
      
    }
    if(signal_i_use==4){
      for(i in 1:length(sequse)){
        entry_function_array[,,i]=hadamard.prod(BBprod[[5]][1,,],X_datacross_reg_SNR[[2]][,,i])+hadamard.prod(BBprod[[6]][1,,],X_datacross_reg_SNR[[3]][,,i])
        entry_function_array[,,i]=entry_function_array[,,i]-hadamard.prod(BBprod[[5]][1,,],matrix(localmean[[2]],64,64))-hadamard.prod(BBprod[[6]][1,,],matrix(localmean[[3]],64,64))
      }
      entry_function_matrix_list[[signal_i_use]]=entry_function_array
    }
    
  }
  
  return(entry_function_matrix_list)
}

true_entry_mean=function(){
  accuracy=0.0000001
  source('./SimDataGeneration/true_signal_different.R')
  #noise level
  sigma=0.1
  adjust_nonlinear=0.3
  #presentation
  adjust_nonlinear=adjust_nonlinear*2
  
  
  sigma_use=c()
  
  sequse=seq(0,1,accuracy)
  
  X_data_SNR=array(0,c(1,1,length(sequse)))
  for(i in 1:1){
    for(j in 1:1){
      X_data_SNR[i,j,]=sequse
    }
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
  
  true_entry_mean_matrix_all=list()
  
  true_entry_mean_matrix=matrix(0,64,64)
  
  
  for(signal_i in 1:4){
    print(signal_i)
    for(i in 1:64){
      print(i)
      for(j in 1:64){
        if(signal_i==1){
          true_entry_mean_matrix[i,j]=mean(BBprod[[1]][1,i,j]*X_datacross_reg_SNR[[1]])
        }
        if(signal_i==2){
          true_entry_mean_matrix[i,j]=mean(BBprod[[2]][1,i,j]*X_datacross_reg_SNR[[2]])
        }
        if(signal_i==3){
          true_entry_mean_matrix[i,j]=mean(BBprod[[3]][1,i,j]*X_datacross_reg_SNR[[2]]+BBprod[[4]][1,i,j]*X_datacross_reg_SNR[[2]])
          
        }
        if(signal_i==4){
          true_entry_mean_matrix[i,j]=mean(BBprod[[5]][1,i,j]*X_datacross_reg_SNR[[2]]+BBprod[[6]][1,i,j]*X_datacross_reg_SNR[[3]])
        }
      }
    }
    true_entry_mean_matrix_all[[signal_i]]=true_entry_mean_matrix
  }
  
  return(true_entry_mean_matrix_all)
  
}

estimated_entryfunction=function(accuracy=0.001,result, SNR=0.1){
  #accuracy=0.00010001

  sigma_use=c()
  
  sequse=seq(0,1,accuracy)
  X_data_SNR=array(0,c(64,64,length(sequse)))
  for(i in 1:64){
    for(j in 1:64){
      X_data_SNR[i,j,]=sequse
    }
  }
  
  
  
  entry_function_matrix_list_iter_list=list()
  entry_function_matrix_list=list()
  entry_function_array=array(0,c(64,64,length(sequse)))
  
  
  #X_truncated_knots_vali=BX_sample_fun(X_data_SNR,order=4,K=6)
  #X_truncated_vali=X_truncated_knots_vali$BX_sample_tensor
  #knots=X_truncated_knots_vali$knots
  #rm(X_truncated_knots_vali)
  num_knots=6
  knots= quantile(c(X_data_SNR), probs = c(seq(0, 1, 1/(num_knots - 1))))
  X_truncated_vali=tildePhiX_trans(X_data_SNR, knots)
  
  middle_list=array(0,c(64,64,6))
  
  for(iter in 1:50){
    print('iter')
    print(iter)
    for(signal_i in 1:4){
      print(signal_i)
      
      BB=result[[iter]][[signal_i]]
      b0=result[[iter]][[signal_i+4]]
      for(i in 1:length(sequse)){
        for(middle_in in 1:6){
          middle_list[,,middle_in]=hadamard.prod(BB[,,middle_in],X_truncated_vali[,,middle_in,i])
        }
        entry_function_array[,,i]=apply(middle_list,c(1,2),FUN=sum)
      }
      entry_function_matrix_list[[signal_i]]=entry_function_array
    }
    entry_function_matrix_list_iter_list[[iter]]=entry_function_matrix_list
    
  }
  
  return(entry_function_matrix_list_iter_list)
}


linear_estimated_entryfunction=function(accuracy=0.001,TLR=1,result=NA,result_vec=NA,id=NA){
  sequse=seq(0,1,accuracy)
  
  X_data_SNR=array(0,c(64,64,length(sequse)))
  for(i in 1:64){
    for(j in 1:64){
      X_data_SNR[i,j,]=sequse
    }
  }
  entry_function_matrix_list=list()
  entry_function_matrix_list_iter_list=list()
  entry_function_array=array(0,c(64,64,length(sequse)))
  
  
  if(TLR==1){
    
    for(iter in 1:50){
      print('iter')
      print(iter)
      for(signal_i in 1:4){
        print(signal_i)
        
        BB=result[[iter]][[signal_i]]
        b0=result[[iter]][[signal_i+4]]
        
        
        for(i in 1:length(sequse)){
          entry_function_array[,,i]=hadamard.prod(BB,X_data_SNR[,,i])
        }
        
        entry_function_matrix_list[[signal_i]]=entry_function_array
      }
      entry_function_matrix_list_iter_list[[iter]]=entry_function_matrix_list
      
    }
    
    return(entry_function_matrix_list_iter_list)
    
  }else{
    
    
    for(iter in 1:50){
      print('Enet')
      print(iter)
      for(signal_i in 1:4){
        print(signal_i)
        
        BB_vec=result_vec[[iter]][[signal_i]]
        b0=result_vec[[iter]][[signal_i+4]]
        BB=matrix(BB_vec,64,64)
        
        for(i in 1:length(sequse)){
          entry_function_array[,,i]=hadamard.prod(BB,X_data_SNR[,,i])
        }
        
        entry_function_matrix_list[[signal_i]]=entry_function_array
      }
      entry_function_matrix_list_iter_list[[iter]]=entry_function_matrix_list
      
    }
    
    return(entry_function_matrix_list_iter_list)
  }
}

linear_integration_of_entry_fun=function(TLR=1,BB=NA,BB_vec=NA){
  if(TLR==1){
    BB=0.5*BB
    return(BB)
  }else{
    BB=matrix(BB_vec,64,64)  
    BB=0.5*BB
    return(BB)
  }
}


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

