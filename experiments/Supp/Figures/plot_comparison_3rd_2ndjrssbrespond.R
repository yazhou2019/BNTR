#Plot for comparison
library(BNTR)
##All the plot
p=c(64,64)
par(mar = c(0,0,0,0),mfcol=c(5,5),mai=c(0.05,0.2,0.3,0.01))

#True signal 
source('./SimDataGeneration/true_signal_different.R')
source('~/Desktop/Github_test/BroadcasTR/clean_code/FullSimPaper/SimDataGeneration/true_signal_different_butterfly.R', echo=TRUE)
index = 2

color_used = c(5/255,120/255,5/255)

BBplot_true = array(1, c(64,64,3))
BBplot = array(1, c(64,64,3))
for(signal_i in 1:length(BB)){
  title=str_c('Case ',signal_i)
  plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
  if(length(dim(BBplot_true)) ==2){
    BBplot_true = 1 - BB[[signal_i]]
  }else{
  BBplot_true[,,-index] = 1 - BB[[signal_i]]
  }
  BBplot_true[,,1] = color_used[1]  * ((-1/color_used[1] + 1) *BB[[signal_i]] + 1/color_used[1])
  BBplot_true[,,2] = color_used[2]  * ((-1/color_used[2] + 1) *BB[[signal_i]] + 1/color_used[2]) 
  BBplot_true[,,3] = color_used[3] * ((-1/color_used[3] + 1) *BB[[signal_i]] + 1/color_used[3])
  rasterImage(BBplot_true, 1, 1, p[1], p[2])
  signal=c('True Signal')
  if(signal_i==1){
    mtext(signal,3,line=0.2)
  }
  mtext(title,2,line=0.2) 
}


method_all=c("TLR1","TLR2","ENet","BNTR")

largest_value = matrix(0, 5, 4)



res_mat_tem = matrix(0, 50, 5)
result_tem = list()
n=1000
for(method_i in 1:4){
  method=method_all[method_i]

  
  if(method_i <=3){
    source('./CollectISEs/ISE_comput_functions.R')
    path1=str_c("./SimResults/", method, "_", n, "_new.Rdata")
    path2=str_c("./SimResults/final_error_", method, "_", n, "_new.Rdata")
    load(path1)
    load(path2)
    
    
    res_mat_tem[,1:4] = t(final_error)
    for(iter in 1:50){
      temtem = list()
      for(signal_i in 1:4){
        temtem[[signal_i]] = result[[iter]][[signal_i]]
      }
      result_tem[[iter]] = temtem 
    }
    
  }else{
    load("./SimResults/BNTR1000_new_K7_2rd.Rdata")  
    load("./SimResults/res_mat1000K8.Rdata")
    
    
    res_mat_tem[,1:4] = res_mat
    for(iter in 1:50){
      temtem = list()
      for(signal_i in 1:4){
        temtem[[signal_i]] = result[[iter]][[signal_i]]
      }
      result_tem[[iter]] = temtem 
    }
    
  }
  
  
  
  path1=str_c("./SimResults/", method, "_", n, "_butterfly.Rdata")
  path2=str_c("./SimResults/res_mat", n, "_", method, "_butterfly.Rdata")
  load(path1)
  load(path2)
  for(iter in 1:50){
    result_tem[[iter]][[5]] = result[[iter]][[1]]
  }
  res_mat_tem[,5] = res_mat[,5]
  
  res_mat = res_mat_tem
  result  = result_tem  
  
  #load the data
  if(method=="TLR1"){
    signal="TLR"
  }
  if(method=="TLR2"){
    signal="TLR-rescaled"
  }
  if(method=="ENet"){
    signal="ENetR"
  }
  if(method=="BNTR"){
    #knots=c(0,0.25,0.5,0.75,1)
    knots = c(0,0.2,0.4,0.6,0.8,1) 
    signal="BroadcasTR"
  }
  
  
  
  
  for(signal_i in 1:5){
      iditer=which(order(res_mat[,signal_i])==26)
      
   
    for(iter in iditer){
      
      plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
      if(method=="BNTR"){
        BB <- fhatnorm_ten(result[[iter]][[signal_i]], knots)
        }else{
        BB=result[[iter]][[signal_i]]  
      }
      if(method!="BNTR"){
        if(sum(BB)!=0){
          BB=sqrt(BB^2/12)
        }
      }
      
      largest_value[signal_i,method_i]=max(BB)
      
      if(length(dim(BBplot)) ==2){
        BBplot = 1 - BB
      }else{
        BBplot[,,-index] = 1 - BB
      }
      #BBplot[,,index] = BB
     # BBplot[,,2] =   BBplot[,,2] *0.5
      BBplot[,,1] = color_used[1]  * ((-1/color_used[1] + 1) *BB + 1/color_used[1])
      BBplot[,,2] = color_used[2]  * ((-1/color_used[2] + 1) *BB + 1/color_used[2]) 
      BBplot[,,3] = color_used[3] * ((-1/color_used[3] + 1) *BB + 1/color_used[3])
      
      rasterImage(BBplot, 1, 1, p[1], p[2]) 
      if(signal_i==1){
        mtext(signal,3,line=0.2)
      }
      
    }
  }
}

#17.38909,111.01420  ,14.839368
