



obtain_linear_BB <-function(){
  library(stringr)
  methods = c("TLR1", "TLR2", "ENet")
  env_name = "/Users/zhouya/Desktop/Github_test/BroadcasTR/clean_code/FullSimPaper/SimResults/"
  sample_sizes = c(500, 750, 1000)
  sample_size_list = list()
  for(size_i in 1:3){
    sample_size = sample_sizes[size_i]
    method_list = list()
    for(method_i in 1:3){
      signal_list  = list()
      data_name = str_c(methods[method_i], "_", sample_size, "_new.Rdata")
      file_name = str_c(env_name, data_name)
      ise_data_name = str_c("final_error_", data_name)
      ise_file_name = str_c(env_name, ise_data_name)
      load(file_name)
      load(ise_file_name)
      res_mat =  t(final_error)
      for(signal_i in 1:4){
        iditer= which(order(res_mat[,signal_i])==26)
        BB = result[[iditer]][[signal_i]]
        signal_list[[signal_i]] = BB
      }
      data_name_bu = str_c(methods[method_i], "_", sample_size, "_butterfly.Rdata") 
      ise_data_name_bu = str_c("res_mat", sample_size, "_", methods[method_i], "_butterfly.Rdata")
      file_name_bu = str_c(env_name, data_name_bu)
      ise_file_name_bu = str_c(env_name, ise_data_name_bu)
      load(file_name_bu)
      load(ise_file_name_bu)
      iditer= which(order(res_mat[,5])==26)
      BB = result[[iditer]][[1]]
      signal_list[[5]] = BB
      method_list[[method_i]] = signal_list
    }
    sample_size_list[[size_i]] = method_list
  }
  return(sample_size_list)
}


obtain_the_curves_values <- function(n_i=4, signal_i=5, distri="unif",qseq_used = c(0.25, 0.75),index=NA, iditer=NA){
  setwd("~/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults")
  library(msm)
  source('./additional_utility.R', echo=TRUE)
  
  
  n_all=c(250, 500, 750, 1000)
  library(stringr)
  library(BNTR)
  
  
  n=n_all[n_i]
  n_size=str_c("n=",n)
  if(distri =="tnorm"){
    if(n==250){
      path1 = "./SimResults/BNTR250_tnorm_K4.Rdata"
      path2 = "./MSE_ISE/res_mat250K4_tnorm_exp.Rdata"
    }
    if(n==500){
      path1 = "./SimResults/BNTR500_tnorm_K5.Rdata"
      path2 = "./MSE_ISE/res_mat500K5_tnorm_exp.Rdata"
    }
    if(n==750){
      path1 = "./SimResults/BNTR750_tnorm_K5.Rdata"
      path2 = "./MSE_ISE/res_mat750K5_tnorm_exp.Rdata"
    }
    if(n==1000){
      path1 = "./SimResults/BNTR1000_tnorm_K7.Rdata"
      path2 = "./MSE_ISE/res_mat1000K7_tnorm_exp.Rdata"
    }
    load(path1)
    load(path2)
  }
  
  
  if(distri =="unif"){
    method = "BNTR"
    if(signal_i <= 4){
      if(n <1000&n>250){
        path1=str_c("/Users/zhouya/Desktop/Github_test/BroadcasTR/clean_code/FullSimPaper/SimResults/", method, "_", n, "_new.Rdata")
        path2=str_c("/Users/zhouya/Desktop/Github_test/BroadcasTR/clean_code/FullSimPaper/SimResults/final_error_", method, "_", n, "_new.Rdata")
        load(path1)
        load(path2)
        res_mat = t(final_error)
      }else if(n ==1000){
        load("/Users/zhouya/Desktop/Github_test/BroadcasTR/clean_code/FullSimPaper/SimResults/BNTR1000_new_K7_2rd.Rdata")  
        load("/Users/zhouya/Desktop/Github_test/BroadcasTR/clean_code/FullSimPaper/SimResults/res_mat1000K8.Rdata")
      }else if(n==250){
        load("/Users/zhouya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/SimResults/BNTR250_new_K4.Rdata")
        load("/Users/zhouya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/MSE_ISE/res_mat250K4uniform.Rdata")
      }
    }
    if(signal_i ==5){
      if(n==1000){
        load("./SimResults/BNTR1000_butterfly_K7.Rdata")
        load("./MSE_ISE/res_mat1000K7butterfly.Rdata")
      }
      if(n==750){
        load("./SimResults/BNTR750_butterfly_K5.Rdata")
        load("./MSE_ISE/res_mat750K5butterfly.Rdata")
      }
      if(n==500){
        load("./SimResults/BNTR500_butterfly_K5.Rdata")
        load("./MSE_ISE/res_mat500K5butterfly.Rdata")
      }
      if(n==250){
        load("./SimResults/BNTR250_butterfly_K4.Rdata")
        load("./MSE_ISE/res_mat250K4uniform.Rdata")
      }
      
    }
    
  }
  if(is.na(iditer)==1){
  iditer= which(order(res_mat[,signal_i])==26)
  }
  
  
  if(distri=="unif" & signal_i ==5){
    BBB=result[[iditer]][[1]]
  }else{
    BBB=result[[iditer]][[signal_i]]
  }
  
  
  tem_tem = plot_fhatminusftrue_vary_sample(signal_i=signal_i, iter=iditer, 
                                        qseq = qseq_used, result=NA, knots = NA, 
                                        n_data=n, onlyone=1,BBB=BBB,distri=distri,index=index)
  tem = tem_tem$res
  x_seq = tem[[1]]$data$iteration[1:100]
  y_true_seq = tem[[1]]$data$ob_value[1:100]
  y_est = list()
  for(res_i in 1:length(tem)){
    y_est[[res_i]] = tem[[res_i]]$data$ob_value[101:200]
  }
  
  final_res = list(x=x_seq, y_true=y_true_seq, y_est = y_est, iditer=iditer, index=tem_tem$index)
  return(final_res)
}


obtain_the_curves_values_linear <- function(size_i=1, signal_i=3, method_i=2, qseq_used=c(0.5), index=NA){
  
  sample_size_list = obtain_linear_BB()
  
  tem_index =   plot_fhatminusftrue_vary_sample_linear(
    sample_size_list,
    signal_i=signal_i, 
    qseq = qseq_used, 
    size_i=size_i, 
    method_i=method_i,
    onlyone=1,
    distri="unif",
    index = index)
  
  tem=tem_index$res
  x_seq = tem[[1]]$data$iteration[1:100]
  y_true_seq = tem[[1]]$data$ob_value[1:100]
  y_est = list()
  for(res_i in 1:length(tem)){
    y_est[[res_i]] = tem[[res_i]]$data$ob_value[101:200]
  }
  
  final_res = list(x=x_seq, y_true=y_true_seq, y_est = y_est, index =  tem_index$index)
  return(final_res)
  
}



#res = obtain_the_curves_values(n_i=4, signal_i=5, distri="unif",qseq_used = c(0.5),index=NA)



obtain_curve_values_all <- function(qseq_used=c(0)){
  
  size_list = list()
  for(size_i in 1:3){
    signal_list = list()
    for(signal_i in 1:5){
      method_list = list()
      use_index = NA
      for(method_i in c(2,1,3)){
        if(method_i <=3){
          res_tem = obtain_the_curves_values_linear(size_i=size_i, signal_i=signal_i, method_i=method_i, qseq_used=qseq_used, index = use_index)  
          if(method_i ==2){
            use_index = res_tem$index
            use_x = res_tem$x
            use_y_true = res_tem$y_true
          }
        }
        method_list[[method_i]] = res_tem
      }
      res_tem = obtain_the_curves_values(n_i=size_i+1, signal_i=signal_i, qseq_used=qseq_used,index=use_index)
      method_list[[4]] = res_tem
      
      signal_list[[signal_i]] = method_list
    }
    
    size_list[[size_i]] = signal_list
  }
  return(size_list)
  
}



obtain_curve_values_BNTR <- function(qseq_used=c(0.5)){
  
  size_list = list()
  for(size_i in 1:3){
    signal_list = list()
    for(signal_i in 1:5){
      method_list = list()
      
      use_index = NA
      res_tem = obtain_the_curves_values(n_i=size_i+1, signal_i=signal_i, qseq_used=qseq_used,index=use_index,iditer = NA)
      use_iter = res_tem$iditer
      use_index = res_tem$index
      method_list[[use_iter]] = res_tem
      
      for(iter_i in c(1:50)[-use_iter]){
      res_tem = obtain_the_curves_values(n_i=size_i+1, signal_i=signal_i, qseq_used=qseq_used,index=use_index,iditer = iter_i)
      method_list[[iter_i]] = res_tem
      }
      
      signal_list[[signal_i]] = method_list
    }
    
    size_list[[size_i]] = signal_list
  }
  return(size_list)
  
}

obtain_curve_values_BNTR_one <- function(qseq_used=c(0.5)){
  
  size_list = list()
  for(size_i in 1:3){
    signal_list = list()
    for(signal_i in 1:5){
      method_list = list()
      
      use_index = NA
      res_tem = obtain_the_curves_values(n_i=size_i+1, signal_i=signal_i, qseq_used=qseq_used,index=use_index,iditer = NA)
      use_iter = res_tem$iditer
      use_index = res_tem$index
      method_list[[1]] = res_tem
      
      signal_list[[signal_i]] = method_list
    }
    
    size_list[[size_i]] = signal_list
  }
  return(size_list)
  
}

size_list = obtain_curve_values_BNTR_one(qseq_used=c(0.5))



obtain_ggplot <- function(size_i=3, size_list){
  method_names = c("TLR", "TLR-rescaled", "ENetR", "BroadcasTR")
  size_names = c("n=500", "n=750", "n=1000")
  small_index_tem = c(seq(1,100,5),100)
  small_index = small_index_tem
  for(i in 2:5){
    small_index = c(small_index, small_index_tem+100*(i-1))
  }
  
  res_signal = list()
  
  for(signal_i in 1:5){
    ob_value = size_list[[1]][[signal_i]][[2]]$y_true
    iteration = c(rep(size_list[[1]][[signal_i]][[2]]$x, 5))
    num_xpoints = length(size_list[[1]][[signal_i]][[2]]$x)
    group = c(rep('True', num_xpoints))
    
    for(method_i in 1:length(method_names)){
      ob_value = c(ob_value, size_list[[size_i]][[signal_i]][[method_i]]$y_est[[1]])
      group = c(group, rep(method_names[method_i], num_xpoints))
    }  
    
    df_small = data.frame(iteration=iteration[small_index] , ob_value=ob_value[small_index], group=group[small_index])
    df = df_small
    res = ggplot(df,aes(x=iteration, y=ob_value,colour=group, linetype=group, group=group))+geom_line(lwd = 0.5) + labs(x=NULL, y="") #+labs(x=NULL,y=str_c("Case ", signal_i)) 
    res = res + scale_color_manual(name="Method",values=c("red","blueviolet","green","orange","black"))
    # res = res + scale_shape_manual("", values=c(4,5,6,10,7))
    res = res + theme(legend.title = element_text(size=16),legend.text = element_text(size=12))
    res = res + scale_linetype_manual(name="Method",values = c("longdash","dotted", "dotdash", "dashed","solid")) 
    res = res + ylim(min(ob_value)-0.1,max(ob_value)+0.1) #+geom_point() #+ labs(title = size_names[size_i]) + theme(plot.title=element_text(hjust=0.5))
    res_signal[[signal_i]] = res
  }
  return(res_signal)
}


obtain_ggplot_BNTR <- function(size_i=3, size_list){
  method_names = c("BroadcasTR")
  size_names = c("n=500", "n=750", "n=1000")
  small_index_tem = c(seq(1,100,5),100)
  small_index = small_index_tem
  #for(i in 2){
  #  small_index = c(small_index, small_index_tem+100*(i-1))
  #}
  
  res_signal = list()
  
  for(signal_i in 1:5){
    ob_value = size_list[[1]][[signal_i]][[2]]$y_true
    iteration = c(rep(size_list[[1]][[signal_i]][[2]]$x, 1))
    num_xpoints = length(size_list[[1]][[signal_i]][[2]]$x)
    group = c(rep('True', num_xpoints))
    df_matrix = matrix(0, length(small_index),55)
    df_matrix[,1] = iteration[small_index]
    df_matrix[,2] = ob_value[small_index]
    df_list = list()
    df_list[[1]]=iteration[small_index]
    df_list[[2]] = ob_value[small_index]
    for(iter_i in 1:50){
      #if(iter_i ==1){
      #ob_value = c(ob_value, size_list[[size_i]][[signal_i]][[iter_i]]$y_est[[1]])
      #group = c(group, rep(method_names[1], num_xpoints))
      #df_list[[1]] = ob_value[small_index]
      #df_list[[2]] = iteration[small_index_tem]
      #df_list[[3]] = c(iteration, iteration)[small_index]
      #}else{
      #df_list[[iter_i + 2]] = size_list[[size_i]][[signal_i]][[iter_i]]$y_est[[1]][small_index_tem]
      #}
      df_list[[iter_i + 2]] = size_list[[size_i]][[signal_i]][[iter_i]]$y_est[[1]][small_index]
      df_matrix[, iter_i + 2] = size_list[[size_i]][[signal_i]][[iter_i]]$y_est[[1]][small_index]
    }  
    df_matrix[,54] = 1
    df_matrix[,55] = 2
    
    df_small = data.frame(df_matrix)
    df_small$X54=rep("True",1)
    df_small$X55=rep("BroadcasTR",1)
   # df_small = data.frame(iteration=iteration[small_index] , ob_value=ob_value[small_index], group=group[small_index])
 

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


obtain_ggplot_BNTR_one <- function(size_i=3, size_list){
  method_names = c("BroadcasTR")
  size_names = c("n=500", "n=750", "n=1000")
  small_index_tem = c(seq(1,100,5),100)
  small_index = small_index_tem
  #for(i in 2){
  #  small_index = c(small_index, small_index_tem+100*(i-1))
  #}
  
  res_signal = list()
  
  for(signal_i in 1:5){
    ob_value = size_list[[1]][[signal_i]][[1]]$y_true
    iteration = c(rep(size_list[[1]][[signal_i]][[1]]$x, 1))
    num_xpoints = length(size_list[[1]][[signal_i]][[1]]$x)
    group = c(rep('True', num_xpoints))
    df_matrix = matrix(0, length(small_index),55)
    df_matrix[,1] = iteration[small_index]
    df_matrix[,2] = ob_value[small_index]
    df_list = list()
    df_list[[1]]=iteration[small_index]
    df_list[[2]] = ob_value[small_index]
    for(iter_i in 1:50){
      #if(iter_i ==1){
      #ob_value = c(ob_value, size_list[[size_i]][[signal_i]][[iter_i]]$y_est[[1]])
      #group = c(group, rep(method_names[1], num_xpoints))
      #df_list[[1]] = ob_value[small_index]
      #df_list[[2]] = iteration[small_index_tem]
      #df_list[[3]] = c(iteration, iteration)[small_index]
      #}else{
      #df_list[[iter_i + 2]] = size_list[[size_i]][[signal_i]][[iter_i]]$y_est[[1]][small_index_tem]
      #}
      df_list[[iter_i + 2]] = size_list[[size_i]][[signal_i]][[1]]$y_est[[1]][small_index]
      df_matrix[, iter_i + 2] = size_list[[size_i]][[signal_i]][[1]]$y_est[[1]][small_index]
    }  
    df_matrix[,54] = 1
    df_matrix[,55] = 2
    
    df_small = data.frame(df_matrix)
    df_small$X54=rep("True",1)
    df_small$X55=rep("BroadcasTR",1)
    # df_small = data.frame(iteration=iteration[small_index] , ob_value=ob_value[small_index], group=group[small_index])
    
    
    df = df_small
    #res = ggplot(df,aes(x=iteration, y=ob_value,colour=group, linetype=group, group=group))+geom_line(lwd = 0.5) + labs(x=NULL, y="") #+labs(x=NULL,y=str_c("Case ", signal_i)) 
    alpha_value = 1
    #res = ggplot(df,aes(x=eval(parse(text='X1')), y=eval(parse(text='X2')),group=eval(parse(text='X54')),linetype=eval(parse(text='X54')),color=eval(parse(text='X54'))))+geom_line(lwd = 0.5) + labs(x=NULL, y="") 
    res = ggplot(df,aes(x=X1, y=X2,group=factor(X54),linetype=X54,color=X54))+ labs(x=NULL, y="") #+geom_line(lwd = 0.5) 
    
    res = res + geom_line(data=df,aes(y=X3,group=X55,linetype=X55,color=X55),alpha=alpha_value)
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


res_size = list()
library(ggplot2)
#size_list = obtain_curve_values_all(qseq_used=c(0))
for(size_i in 1:3){
  res_size[[size_i]] = obtain_ggplot_BNTR(size_i=size_i, size_list) 
}

res_size = list()
library(ggplot2)
#size_list = obtain_curve_values_all(qseq_used=c(0))
for(size_i in 1:3){
  res_size[[size_i]] = obtain_ggplot_BNTR_one(size_i=size_i, size_list) 
}


library(ggpubr)
library(stringi)
library(stringr)
png("temtem.png",units="in", width=9.52, height=11,res=300)
ggarrange(res_size[[1]][[1]]+labs(x=NULL,y=str_c("Case ", 1))+ 
            labs(title = "n=500") + theme(plot.title=element_text(hjust=0.5,margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
          res_size[[2]][[1]]+ labs(title = "n=750") + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
          res_size[[3]][[1]]+ labs(title = "n=1000") + theme(plot.title=element_text(hjust=0.5,margin=margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))),
          res_size[[1]][[2]]+labs(x=NULL,y=str_c("Case ", 2)),res_size[[2]][[2]],res_size[[3]][[2]],
          res_size[[1]][[3]]+labs(x=NULL,y=str_c("Case ", 3)),res_size[[2]][[3]],res_size[[3]][[3]],
          res_size[[1]][[4]]+labs(x=NULL,y=str_c("Case ", 4)),res_size[[2]][[4]],res_size[[3]][[4]],
          res_size[[1]][[5]]+labs(x=NULL,y=str_c("Case ", 5)),res_size[[2]][[5]],res_size[[3]][[5]],
          ncol=3, nrow=5, common.legend = TRUE, legend='right')

dev.off()




