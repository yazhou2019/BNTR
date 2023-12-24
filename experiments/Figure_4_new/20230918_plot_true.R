library(stringi)
library(stringr)

setwd("~/Desktop/Research/JRSSB/换editor上传/Table_1/BroadcasTR")

#setwd("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper")
source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')
source('./SimDataGeneration_He/20230224_get_outputs.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_ise.R', echo=FALSE)
#source("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimDataGeneration_He/20230224_get_BB.R")
source("./SimDataGeneration_He/20230224_get_BB.R")

setwd("~/Desktop/Research/JRSSB/换editor上传/Table_1")
#Plot for comparison
library(BNTR)
##All the plot
p=c(64,64)

png("com_all_new_editor.png",units="in", width=8, height=8,res=300)

par(mar = c(0,0,0,0),mfcol=c(5,5),mai=c(0.05,0.2,0.3,0.01))

plot_scale_par = 2.5


index = 3

color_used = c(5/255,120/255,5/255) 

BBplot_true = array(1, c(64,64,3))
BBplot = array(1, c(64,64,3))
BBprod = get_BB()
BB=get_true_l2norm(BBprod)
for(signal_i in c(1,3,2,4,5)){
  title=str_c('Case ',signal_i)
  plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
  BB_signal_i_normlized = plot_scale_par * BB[[signal_i]] #abs(BB[[signal_i]])/max(abs(BB[[signal_i]]))
  
  BB_signal_i_normlized[BB_signal_i_normlized>1] = 1
  
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
  if(signal_i ==2){
    title = "Case 3"
  }
  if(signal_i==3){
    title="Case 2"
  }
  mtext(title,2,line=0.2) 
}


method_all=c("TLR1", "TLR2","ENet","BNTR")
#method_all=c("BNTR")
largest_value = matrix(0, 5, 4)



res_mat_tem = matrix(0, 50, 5)
result_tem = list()
n=1000

methods = c("TLR1", "TLR-rescaled","ENetR", "BroadcasTR")

for(method_i in 1:4){
  method=method_all[method_i]
  # 
  # 
  # #load the data
  # if(method=="TLR1"){
  #  
  #   # path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/TLR1/ise_TLR1_",n,"_20230302.Rdata")
  #   # path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/TLR1/TLR1_", n, "_20230302.Rdata")
  #   path1 = str_c("./TLR1/TLR1_",n,"_20230302.Rdata")
  #   path2 = str_c("./TLR1/ise_TLR1_", n, "_20230302.Rdata")
  #   
  #   signal="TLR"
  # }
  # 
  # 
  # if(method =="BNTR"){
  #   # path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/BNTR",n,"_new_K8_20230302.Rdata")
  #   # path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_BNTR", n, "_new_K8_20230302.Rdata")
  #   path1 = str_c("./SimResults/BNTR",n,"_new_K8_20230302.Rdata")
  #   path2 = str_c("./SimResults/ise_BNTR", n, "_new_K8_20230302.Rdata")
  #   
  #   signal = "BroadcasTR"
  # }
  # 
  # if(method=="TLR2"){
  #   # path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/TLR2_",n, "_linear_new_20230302.Rdata")
  #   # path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_linear",n,"_new_K8_20230302.Rdata")
  #   path1 = str_c("./SimResults/TLR2_",n, "_linear_new_20230302.Rdata")
  #   path2 = str_c("./SimResults/ise_linear",n,"_new_K8_20230302.Rdata")
  #   
  #   signal="TLR-rescaled"
  # }
  # 
  # if(method=="ENet"){
  #   # path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ENet_",n,"_new_20230302.Rdata")
  #   # path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_ENet", n,"_new_K8_20230302.Rdata")
  #   path1 = str_c("./SimResults/ENet_",n,"_new_20230302.Rdata")
  #   path2 = str_c("./SimResults/ise_ENet", n,"_new_K8_20230302.Rdata")
  #   
  #   signal = "ENetR"
  # }
  # 
  # load(path1)
  # load(path2)
  
  
  method_prime = method
  signal = methods[method_i]
  path1 = str_c("./", signal, "/SimResults/ise_",method_prime,"_", n, "_20230918.Rdata" )
  load(path1)
  path2 = str_c("./", signal, "/SimResults/",method_prime,"_", n, "_20230918.Rdata" )
  load(path2)
  
  
  
  
  
  
  for(signal_i in c(1,3,2,4,5)){
    iditer= order(ise_mat[,2*signal_i-1])[26]#which(order(ise_mat[,2*signal_i-1])==50)
    #iditer = c(6)
    
    for(iter in iditer){
      
      plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
      
      
      if(method=="BNTR"){
        knots_used = result[[iter]][[signal_i+15]]
        BB = plot_scale_par*fhatnorm_ten(result[[iter]][[signal_i]],knots_used)
        # normalization
        print(c(max(BB),min(BB)))
        #BB = BB/max(BB)
      } else {
        BB=result[[iter]][[signal_i]] 
        if(sum(BB)!=0){
          BB=plot_scale_par* sqrt(BB^2/12)
        }
        # normalization
        print(c(max(BB),min(BB)))
       # BB = BB/max(BB)
      }
      BB[BB>1]=1
      
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

dev.off()
