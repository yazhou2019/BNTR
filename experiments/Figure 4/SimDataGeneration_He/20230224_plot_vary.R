#plot vary size


setwd("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper")
source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')
#source('./SimDataGeneration_new/20221018_get_outputs.R', echo=FALSE)
#source('./SimDataGeneration_new/20221018_get_ise.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_outputs.R', echo=FALSE)
source('./SimDataGeneration_He/20230224_get_ise.R', echo=FALSE)
plot_scale_par = 2.5

library(stringi)
library(stringr)
#plot vary size



##plot for BNTR in various  sample size
p=c(64,64)

png("temtemtem.png",units="in", width=5.4, height=8,res=300)
par(mar = c(0,0,0,0),mfcol=c(5,3),mai=c(0.05,0.25,0.3,0.2))


n_all=c(250, 500, 750, 1000)
library(stringr)
library(BNTR)

BB_mat <- array(0, c(64,64,5,4))
largest_value = matrix(0, 5, 4)

color_used = c(5/255,120/255,5/255)
scale_par = 1


for(n_i in 2:4){
  method="BNTR" # nonlinear model
  #method="TLR2" # linear model with rescaling
  
  #distri="tnorm"
  distri="beta"
  
  n=n_all[n_i]
  n_size=str_c("n=",n)
  
  
  
  
  
  #load the data
  
  path1 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/", method,n,"_new_K8_20230302.Rdata")
  path2 = str_c("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/FullSimPaper/SimResults/ise_", method, n, "_new_K8_20230302.Rdata")
  
  
  load(path1)
  load(path2)
  

  for(signal_i in c(1,3,2,4,5)){
    iditer=order(ise_mat[,2*signal_i-1])[26]
    
    
    for(iter in iditer){
      
      plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
      
      if(method=="BNTR"){
        #BB=BB_l2norm_squre(result[[iter]][[signal_i]],knots,dimension = 2,original=1)-integration_of_entry_fun(result[[iter]][[signal_i]],knots)^2
        #BB= fhatnorm_ten(result[[iter]][[signal_i]],knots)
        
        knots_used = result[[iter]][[signal_i+15]]
        BB = plot_scale_par*fhatnorm_ten(result[[iter]][[signal_i]],knots_used)
        BB[BB>1] = 1
        
        
        BB_mat[,, signal_i, n_i]=BB
      }else {
        BB=result[[iter]][[signal_i]] 
        if(sum(BB)!=0){
          BB=sqrt(BB^2/12)
        }
      }
      
      largest_value[signal_i, n_i]=max(BB)
      
      BBplot_true = array(1, c(64,64,3))
      #BBplot_true[,,-2] = 1 - BB
      
      BBplot_true[,,1] = color_used[1]  * ((-1/color_used[1] + 1) *BB/scale_par + 1/color_used[1])
      BBplot_true[,,2] = color_used[2]  * ((-1/color_used[2] + 1) *BB/scale_par+ 1/color_used[2]) 
      BBplot_true[,,3] = color_used[3] * ((-1/color_used[3] + 1) *BB/scale_par + 1/color_used[3])
      
      #rasterImage(BB, 1, 1, p[1], p[2]) 
      rasterImage(BBplot_true, 1, 1, p[1], p[2]) 
      
      if(signal_i==1){
        mtext(n_size,3,line=0.2)
      }
      if(n==250){
        title=str_c('Case ',signal_i)
        mtext(title,2,line=0.2) 
      }
    }
  }
}
dev.off()

