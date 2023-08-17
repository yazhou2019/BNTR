
source('./CollectISEs/ISE_comput_functions.R')
source('./SimDataGeneration/true_signal_different.R')
library(rTensor)
library(matrixcalc)


n_all=c(500,750,1000)
method_all=c("TLR1", "TLR2", "ENet")

n=n_all[3] # user define
method=method_all[2] #user define



#path=str_c("/Users/zhouya/Desktop/simulation for paper/Simulation_paper_new/results/",method,"_",n,"_new.Rdata")
path=str_c("./SimResults/", method, "_", n, "_linear_new.Rdata")
load(path)

source('./CollectISEs/ISE_comput_functions.R')

entry_function_matrix_list_iter_list=linear_estimated_entryfunction(accuracy=0.001,TLR=1,result=result)

cons_mean_esimation=matrix(0,4,50)

#constant part after de mean of the estimation
for(signal_i in 1:4){
  for(iter in 1:50){
    BB=result[[iter]][[signal_i]]
    b0=result[[iter]][[signal_i+4]]
    cons_mean_esimation[signal_i,iter]=sum(linear_integration_of_entry_fun(TLR=1,BB=BB))+b0
  }
}



true_mean=sum_mean_true_signal(accuracy=0.00010001)
cons_err=matrix(0,4,50)
for(signal_i in 1:4){
  for(iter in 1:50){
    cons_err[signal_i,iter]=cons_mean_esimation[signal_i,iter]-true_mean$summean[[signal_i]]-1
  }
}



#de mean of the estimated function
for(iter in 1:50){
  for(signal_i in 1:4){
    BB=matrix(result[[iter]][[signal_i]],64,64)
    middle=linear_integration_of_entry_fun(TLR=1,BB)
    for(i in 1:1001){
      entry_function_matrix_list_iter_list[[iter]][[signal_i]][,,i]=entry_function_matrix_list_iter_list[[iter]][[signal_i]][,,i]-middle
    }
    
  }
}







#true local sum
true_sum_localmean=sum_mean_true_signal(accuracy=0.00010001)

#has de mean
true_entryfunction_use=true_entryfunction(accuracy=0.001,true_sum_localmean$localmean)

functionpart=matrix(0,4,50)

for(iter in 1:50){
  for(signal_i in 1:4){
    functionpart[signal_i,iter]=sum(1/1001*(entry_function_matrix_list_iter_list[[iter]][[signal_i]]-true_entryfunction_use[[signal_i]])^2)
    
  }
}  



final_error_2=list()
final_error_2[[1]]=cons_err
final_error_2[[2]]=functionpart


final_error=cons_err^2+functionpart



hehe_matrix=matrix(0,2,4)

for(signal_i in 1:4){
  hehe_matrix[1,signal_i]=mean(final_error[signal_i,])
  hehe_matrix[2,signal_i]=sqrt(var(final_error[signal_i,]))
}
round(hehe_matrix,4)
setwd("./SimResults")

name=str_c("final_error_", method, "_", n, "_new.Rdata")
save(final_error, file=name)
setwd("../")



