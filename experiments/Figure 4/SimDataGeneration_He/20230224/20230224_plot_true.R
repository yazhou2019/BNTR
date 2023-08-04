setwd("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper")
source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')
#source('./SimDataGeneration_new/20221018_get_outputs.R', echo=FALSE)
#source('./SimDataGeneration_new/20221018_get_ise.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_outputs.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_ise.R', echo=FALSE)

#source("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/get_true_normtensor.R")

source("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimDataGeneration_He/20230224_get_BB.R")
#Plot for comparison
library(BNTR)
##All the plot
p=c(64,64)
par(mar = c(0,0,0,0),mfcol=c(5,5),mai=c(0.05,0.2,0.3,0.01))




index = 3

color_used = c(5/255,120/255,5/255) 

BBplot_true = array(1, c(64,64,3))
BBplot = array(1, c(64,64,3))
BBprod = get_BB()
BB=get_true_l2norm(BBprod)
for(signal_i in 1:length(BB)){
  title=str_c('Case ',signal_i)
  plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
  BB_signal_i_normlized = 2 * BB[[signal_i]] #abs(BB[[signal_i]])/max(abs(BB[[signal_i]]))
  if(length(dim(BBplot_true)) ==2){
    BBplot_true = 1 - BB_signal_i_normlized
  }else{
    BBplot_true[,,-index] = 1 - BB_signal_i_normlized
  }
  print(c(max(BB[[signal_i]]),min((BB[[signal_i]]))))
  BBplot_true[,,1] = color_used[1]  * ((-1/color_used[1] + 1) *BB_signal_i_normlized + 1/color_used[1])
  BBplot_true[,,2] = color_used[2]  * ((-1/color_used[2] + 1) *BB_signal_i_normlized + 1/color_used[2]) 
  BBplot_true[,,3] = color_used[3] * ((-1/color_used[3] + 1) *BB_signal_i_normlized + 1/color_used[3])
  rasterImage(BBplot_true, 1, 1, p[1], p[2])
  signal=c('True Signal')
  if(signal_i==1){
    mtext(signal,3,line=0.2)
  }
  mtext(title,2,line=0.2) 
}


method_all=c("TLR1", "TLR2","ENet","BNTR")
#method_all=c("BNTR")
largest_value = matrix(0, 5, 4)



res_mat_tem = matrix(0, 50, 5)
result_tem = list()
n=1000

for(method_i in 4){
  method=method_all[method_i]
  
  
  #load the data
  if(method=="TLR1"){
    signal="TLR"
  }
  
  
  if(method =="BNTR"){
   # path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/BNTR",n,"_new_K8_20230103.Rdata")
   # path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_BNTR", n, "_new_K8_20230103_rank1_8.Rdata")
    path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/BNTR",n,"_new_K8_20230224.Rdata")
    path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_BNTR", n, "_new_K8_20230224.Rdata")
    
    signal = "BroadcasTR"
  }
  
  if(method=="TLR2"){
    path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/TLR2_",n, "_linear_new_20230103.Rdata")
    path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_linear",n,"_new_K8_20230103.Rdata")
    signal="TLR-rescaled"
  }
  
  if(method=="ENet"){
    path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ENet_",n,"_new_20230103.Rdata")
    path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_ENet", n,"_new_K8_20221222.Rdata")
    signal = "ENetR"
  }
  
  load(path1)
  load(path2)
  
  
  
  
  
  
  
  for(signal_i in 1:5){
    iditer= which(order(ise_mat[,2*signal_i-1])==49)
    
    
    for(iter in iditer){
      
      plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
      
      
      if(method=="BNTR"){
        knots_used = result[[iter]][[signal_i+15]]
        BB = fhatnorm_ten(result[[iter]][[signal_i]],knots_used)
        # normalization
        print(c(max(BB),min(BB)))
        #BB = BB/max(BB)
      } else {
        BB=result[[iter]][[signal_i]] 
        if(sum(BB)!=0){
          BB=sqrt(BB^2/12)
        }
        # normalization
        print(c(max(BB),min(BB)))
       # BB = BB/max(BB)
      }
      
      #largest_value[signal_i, n_i]=max(BB)
      BBplot = array(1, c(64,64,3))
      BBplot[,,-2] = 1 - BB
      
      BBplot[,,1] = color_used[1]  * ((-1/color_used[1] + 1) *BB + 1/color_used[1])
      BBplot[,,2] = color_used[2]  * ((-1/color_used[2] + 1) *BB + 1/color_used[2]) 
      BBplot[,,3] = color_used[3] * ((-1/color_used[3] + 1) *BB + 1/color_used[3])
      
      rasterImage(BBplot, 1, 1, p[1], p[2]) 
      if(signal_i==1){
        mtext(signal,3,line=0.2)
      }
      #if(n==250){
      #  title=str_c('Case ',signal_i)
      #  mtext(title,2,line=0.2) 
      #}
    }
  }
  
}

#17.38909,111.01420  ,14.839368
