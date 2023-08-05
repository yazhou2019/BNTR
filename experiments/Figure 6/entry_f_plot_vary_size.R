
#setwd("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper")
library(stringi)
library(stringr)
library(ggplot2)

source("./ComponentsNonLin/functions_needed.R")
source("./SimDataGeneration_He/20230224_get_ise.R")
source("./SimDataGeneration_He/20230224_get_BB.R")

obtain_median_BNTR_iter <-function(case = 1, n = 500){
  #path = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/BNTR", n, "_new_K8_20230302.Rdata") 
  #path = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_BNTR",n,"_new_K8_20230302.Rdata")
  path = str_c("./SimResults/ise_BNTR",n,"_new_K8_20230302.Rdata")
  load(path)
  iter_index = order(ise_mat[,2*case-1])[26] #which(order(ise_mat[,2*case-1])==26)
  return(iter_index)
}


obtain_hat_curve <- function(n=500){
  #path = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/BNTR", n, "_new_K8_20230302.Rdata")
  path = str_c("./SimResults/BNTR", n, "_new_K8_20230302.Rdata")
  load(path)
  y_hat_all_list = list()
  for(iter in 1:50){
    y_hat_signal_list = list()
    for(signal_i in 1:5){
      BB = result[[iter]][[signal_i]]
      bb = result[[iter]][[signal_i+5]]
      knots_used = result[[iter]][[signal_i+15]]
      y_hat_signal_list[[signal_i]] = obtain_hat_curve_signal_i(b0=bb, BB=BB, case=signal_i, precision=1000, knots = knots_used, order = 4)
    }
  y_hat_all_list[[iter]] = y_hat_signal_list
  }
  return(y_hat_all_list)
}



obtain_hat_curve_signal_i <- function(b0, BB, case=1, precision=1000, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4){
  fhat_intten = fhat_ten_int(BB,  precision = precision, knots=knots, order = order)
  y_hat_ten = array(0, c(64, 64, 100))
  for(i in 1:100){
    X=array(i/100, c(64,64))
    y_hat_ten[,,i] = fhat_function_ten(X, BB,knots=knots, order =order) - fhat_intten
  }
  return(y_hat_ten)
}

obtain_true_curve <- function(precision=1000){
  load("./SimResults/Simul_n1000_rep50_final_fix_new.Rdata")
  BBv = data_all[[51]]
  vv = 1
  y_true_signal_list = list()
  for(case in 1:5){
  ftrue_ten_int_case_array = get_ftrue_ten_int_case_array(case=case, BBv=BBv, p=c(64,64), precision=precision)
  y_true_ten = array(0, c(64, 64, 100))
  for(i in 1:100){
    X=array(i/100, c(64,64))
    #y_hat_ten[,,i] = fhat_function_ten(X, BB,knots=knots, order =order) - fhat_intten
    y_true_ten[,,i] = ftrue_demean_ten(X, ftrue_ten_int_case_array, BBv, case=case)
  }
  y_true_signal_list[[case]] = y_true_ten
  }
  return(y_true_signal_list)
}


obtain_index_list <- function(y_hat_ten, y_true_ten, n, qtile=0.55555){
  entry_l2_norm_diff = array(0, c(64, 64, 5))
  for(signal_i in 1:5){
  iter_index = obtain_median_BNTR_iter(case = signal_i, n = n)
  y_tem = y_hat_ten[[iter_index]][[signal_i]] - y_true_ten[[signal_i]]
  for(i1 in 1:64){
    for(i2 in 1:64){
      entry_l2_norm_diff[i1,i2,signal_i] = sum(abs(y_tem[i1,i2,])^2 )
    }
  }
  }
  
  BBprod = get_BB()
  BB=get_true_l2norm(BBprod)
  
  index_list = list()
  for(signal_i in 1:5){
  BB_vec = c(BB[[signal_i]][1,,])
  BB_vec_nonzero = BB_vec[BB_vec!=0]
  nonzeroindex = which(BB_vec!=0)
  
  entry_l2norm_vec = c(entry_l2_norm_diff[,,signal_i])
  entry_l2norm_vec_nonzero = entry_l2norm_vec[nonzeroindex]
  
  #temp0 = round(quantile(1:length(nonzeroindex), 0.55555))当前使用的
  temp0 = round(quantile(1:length(nonzeroindex), qtile))
  temp = c()
  for(i in 1:length(temp0)){
    temp[i] = nonzeroindex[order(entry_l2norm_vec_nonzero)[temp0[i]]]
  }
  
  temp2 = array(1:64^2, c(64,64))
    for(i1 in 1:64){
      for(i2 in 1:64){
        if(temp2[i1,i2]==temp[i]){
          index_list[[signal_i]] = c(i1,i2)
          print(c(entry_l2_norm_diff[i1,i2, signal_i], temp0[i]))
        }
      }
    }
    
  }
  
  
  
  
  return(index_list)
  
}

obtain_ggplot_BNTR <- function(used_index_list = c(10,10), y_hat_ten, y_true_ten){
  method_names = c("BroadcasTR")
  size_names = c("n=500", "n=750", "n=1000")
  small_index_tem = 1:100 # c(seq(1,100,5),100)
  small_index = small_index_tem
  res_signal = list()
  
  for(signal_i in 1:5){
    used_index =  used_index_list[[signal_i]]
    ob_value = y_true_ten[[signal_i]][used_index[1], used_index[2], ]
    iteration = 1:100/100
    num_xpoints = 100
    group = c(rep('True', num_xpoints))
    df_matrix = matrix(0, length(small_index),55)
    df_matrix[,1] = iteration[small_index]
    df_matrix[,2] = ob_value[small_index]
    df_list = list()
    df_list[[1]]=iteration[small_index]
    df_list[[2]] = ob_value[small_index]
    for(iter_i in 1:50){
      df_list[[iter_i + 2]] = y_hat_ten[[iter_i]][[signal_i]][used_index[1], used_index[2], ][small_index]
     
      df_matrix[, iter_i + 2] = y_hat_ten[[iter_i]][[signal_i]][used_index[1], used_index[2], ][small_index]
    }  
    df_matrix[,54] = 1
    df_matrix[,55] = 2
    
    df_small = data.frame(df_matrix)
    df_small$X54=rep("True",1)
    df_small$X55=rep("BroadcasTR",1)
  
    
    df = df_small
    #res = ggplot(df,aes(x=iteration, y=ob_value,colour=group, linetype=group, group=group))+geom_line(lwd = 0.5) + labs(x=NULL, y="") #+labs(x=NULL,y=str_c("Case ", signal_i)) 
    alpha_value = 0.2
    #res = ggplot(df,aes(x=eval(parse(text='X1')), y=eval(parse(text='X2')),group=eval(parse(text='X54')),linetype=eval(parse(text='X54')),color=eval(parse(text='X54'))))+geom_line(lwd = 0.5) + labs(x=NULL, y="") 
    res = ggplot(df,aes(x=X1, y=X2,group=factor(X54),linetype=X54,color=X54))+ labs(x=NULL, y="") #+geom_line(lwd = 0.5) 
    
    res = res + geom_line(data=df,aes(y=X3,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X4,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X5,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X6,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X7,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X8,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X9,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X10,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X11,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X12,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X13,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X14,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X15,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X16,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X17,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X18,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X19,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X20,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X21,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X22,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X23,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X24,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X25,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X26,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X27,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X28,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X29,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X30,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X31,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X32,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X33,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X34,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X35,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X36,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X37,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X38,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X39,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X40,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X41,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X42,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X43,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X44,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X45,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X46,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X47,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X48,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X49,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X50,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X51,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X51,group=X55,linetype=X55,color=X55),alpha=alpha_value)
    res = res + geom_line(data=df,aes(y=X2,group=X54,linetype=X54,color=X54),alpha=1,lwd = 0.5)
    #for(iter_i in 1:50){
    #res = res + geom_line(data=df,aes(y=eval(parse(text=str_c('X',iter_i+2))),group=eval(parse(text='X55')),linetype=eval(parse(text='X55')),color=eval(parse(text='X55'))),alpha=alpha_value)   
    #}    
    
    
    ### "#01A001"
    res = res + scale_color_manual(name="Method",values=c("red","black"))
    
    res = res + theme(legend.title = element_text(size=16),legend.text = element_text(size=12))
    res = res + scale_linetype_manual(name="Method",values = c("dashed","solid")) 
    res = res + ylim(min(ob_value)-0.3,max(ob_value)+0.3) #+geom_point() #+ labs(title = size_names[size_i]) + theme(plot.title=element_text(hjust=0.5))
    res_signal[[signal_i]] = res
  }
  return(res_signal)
}



#y_hat_ten = obtain_hat_curve(n=1000)
#y_true_ten = obtain_true_curve()
#index_list = obtain_index_list(y_hat_ten, y_true_ten, n=1000)
#ggres = obtain_ggplot_BNTR(used_index_list =index_list, y_hat_ten, y_true_ten)



if(1==1){
# for adjust ment
y_hat_ten_list = list()

y_true_ten = obtain_true_curve()
n_size = c(500,750,1000)
for(size_i in 1:3){
  y_hat_ten_list[[size_i]] = obtain_hat_curve(n=n_size[size_i])
}





qtile=0.57
index_list= obtain_index_list(y_hat_ten_list[[3]], y_true_ten, n=1000, qtile=qtile)
index_list_new = list()
for(i in c(1,3,4)){
  index_list_new[[i]] = index_list[[i]]
}

for(i in c(2,5)){
  index_list_new[[i]] = index_list[[i]]
}
index_list = index_list_new

res_size = list()
for(size_i in 1:3){
res_size[[size_i]] = obtain_ggplot_BNTR(used_index_list =index_list, y_hat_ten_list[[size_i]], y_true_ten)
}

library(ggpubr)
library(stringi)
library(stringr)
save_name = str_c("tem", qtile, ".png")
#png(save_name,units="in", width=9.52, height=11,res=300)
ggarrange(res_size[[1]][[1]]+labs(x=NULL,y=str_c("Case ", 1))+ 
            labs(title = "n=500") + theme(plot.title=element_text(hjust=0.5,margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
          res_size[[2]][[1]]+ labs(title = "n=750") + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
          res_size[[3]][[1]]+ labs(title = "n=1000") + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
          res_size[[1]][[3]]+labs(x=NULL,y=str_c("Case ", 2)),res_size[[2]][[3]],res_size[[3]][[3]],
          res_size[[1]][[2]]+labs(x=NULL,y=str_c("Case ", 3)),res_size[[2]][[2]],res_size[[3]][[2]],
          res_size[[1]][[4]]+labs(x=NULL,y=str_c("Case ", 4)),res_size[[2]][[4]],res_size[[3]][[4]],
          res_size[[1]][[5]]+labs(x=NULL,y=str_c("Case ", 5)),res_size[[2]][[5]],res_size[[3]][[5]],
          ncol=3, nrow=5, common.legend = TRUE, legend='right')

#dev.off()

}








