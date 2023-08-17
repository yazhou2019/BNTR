


source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')




library(matrixcalc)
source('./CollectISEs/ISE_comput_functions.R')
source('./SimDataGeneration/true_signal_different.R')

n=500

method="BNTR"
ndata = n
isbutterfly = 0
# load the data set
if(isbutterfly == 0){
    if(ndata==1000){
        loadname = str_c("./SimResults/BNTR", ndata, str_c("_", "new_cv_to_old.Rdata"))
    }else if(ndata==500||ndata==750){
        loadname = str_c("./SimResults/BNTR", ndata, str_c("_", "new_cv_to_old.Rdata"))
    }else if(ndata==250){
        loadname = str_c("./SimResults/BNTR", ndata, str_c("_", "new_cv_to_old.Rdata"))
    }
}else if(isbutterfly == 1){
    if(ndata==1000){
        loadname = str_c("./SimResults/BNTR", ndata, '_butterfly_K7.Rdata')
    }else if(ndata==500||ndata==750){
        loadname = str_c("./SimResults/BNTR", ndata, '_butterfly_K5.Rdata')
    }else if(ndata==250){
        loadname = str_c("./SimResults/BNTR", ndata, '_butterfly_K4.Rdata')
    }
}
load(loadname)



entry_function_matrix_list_iter_list=estimated_entryfunction(accuracy=0.001,result=result)


cons_mean_esimation=matrix(0,4,50)
knots=c(0,0.25,0.5,0.75,1)
for(signal_i in 1:4){
  for(iter in 1:50){
    BB=result[[iter]][[signal_i]]
    b0=result[[iter]][[signal_i+4]]
    cons_mean_esimation[signal_i,iter]=sum(integration_of_entry_fun(BB,knots))+b0
  }
}


true_mean=sum_mean_true_signal(accuracy=0.00010001)
cons_err=matrix(0,4,50)


for(signal_i in 1:4){
  for(iter in 1:50){
    cons_err[signal_i,iter]=cons_mean_esimation[signal_i,iter]-true_mean$summean[[signal_i]]-1
  }
}

#de-mean of the estimated function
for(iter in 1:50){
  for(signal_i in 1:4){
    BB=result[[iter]][[signal_i]]
    middle=integration_of_entry_fun(BB,knots)
    for(i in 1:1001){
      entry_function_matrix_list_iter_list[[iter]][[signal_i]][,,i]=entry_function_matrix_list_iter_list[[iter]][[signal_i]][,,i]-middle
    }
    
  }
}


true_sum_localmean=sum_mean_true_signal(accuracy=0.00010001)
true_entryfunction_use=true_entryfunction(accuracy=0.001,true_sum_localmean$localmean)


functionpart=matrix(0,4,50)
for(iter in 1:50){
  for(signal_i in 1:4){
    functionpart[signal_i,iter]=sum(1/1001*(entry_function_matrix_list_iter_list[[iter]][[signal_i]]-true_entryfunction_use[[signal_i]])^2)
    
  }
}  

final_error=cons_err^2+functionpart



final_error_2=list()
final_error_2[[1]]=cons_err
final_error_2[[2]]=functionpart
final_error=final_error_2[[1]]^2+final_error_2[[2]]



#############see the result
summary_ISEs=matrix(0,2,4)

for(signal_i in 2){
  summary_ISEs[1,signal_i]=mean(final_error[signal_i,])
  summary_ISEs[2,signal_i]=sqrt(var(final_error[signal_i,]))
}
round(summary_ISEs,3)


setwd("./SimResults")

name=str_c("final_error_", method, "_", n, "_new_cv.Rdata")
save(final_error, file=name)
setwd("../")
