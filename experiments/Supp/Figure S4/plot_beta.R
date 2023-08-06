library(stringi)
library(stringr)
#plot vary size



##plot for BNTR in various  sample size
p=c(64,64)

#png("temtemtem.png",units="in", width=5.4, height=8,res=300)
par(mar = c(0,0,0,0),mfcol=c(5,3),mai=c(0.05,0.25,0.3,0.2))


n_all=c(250, 500, 750, 1000)
library(stringr)
library(BNTR)

BB_mat <- array(0, c(64,64,5,4))
largest_value = matrix(0, 5, 4)

color_used = c(5/255,120/255,5/255)
scale_par = 1.1


for(n_i in 2:4){
  method="BNTR" # nonlinear model
  #method="TLR2" # linear model with rescaling
  
  #distri="tnorm"
  distri="beta"
  
  n=n_all[n_i]
  n_size=str_c("n=",n)
  
  
  if(method =="BNTR" && distri=="tnorm"){
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
  }
  
  if(method == "TLR2" && distri=="tnorm"){
    path1 = str_c("./SimResults/TLR2_",n, "_linear_new_tnorm.Rdata")
    path2 = str_c("./MSE_ISE/res_mat", n, "_linear_new_tnorm.Rdata")
  }
  
  
  if(method =="BNTR" && distri=="beta"){
    if(n==500){
      path1 = "./SimResults/BNTR500_beta_K5_2_2.Rdata"
      #path2 = "./MSE_ISE/res_mat500K5_tnorm_exp.Rdata"
    }
    if(n==750){
      path1 = "./SimResults/BNTR750_beta_K5_2_2.Rdata"
      #path2 = "./MSE_ISE/res_mat750K5_tnorm_exp.Rdata"
    }
    if(n==1000){
      path1 = "./SimResults/BNTR1000_beta_K7_2_2.Rdata"
      #path2 = "./MSE_ISE/res_mat1000K7_tnorm_exp.Rdata"
    }
  }
  
  
  
  
  load(path1)
  if(distri!="beta"){
    load(path2)
  }
  if(distri=="beta"){
    
    temfunction <-function(){
      path2 = "./SimResults/ise_meta_beta_n500_750_1000_2_2.Rdata" 
      load(path2)
      return(result)
    }  
    result_tem = temfunction()
    
    res_mat_tem = result_tem[[n_i-1]]
    res_mat = matrix(0, 50, 5)
    for(tem_i in 1:5){
      res_mat[,tem_i] = res_mat_tem[[tem_i]][,3]
    }
  }
  
  
  
  
  
  
  #load the data
  if(method=="BNTR"){
    
    if(n==1000){
      num_knots <- 6
    }else if(n==500||n==750){
      num_knots <- 5
    }else if(n==250){
      num_knots <- 4
    }
    
    if(distri=="tnorm"){
      X_tem = array(rtnorm(64*64*1000, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,1000))
      knots = quantile(c(X_tem), probs = c(seq(0, 1, 1/(num_knots - 1)))) 
    }
    
    if(distri =="uniform"){
      X_tem = array(runif(64*64*1000,0,1),c(64,64,1000))
      knots = quantile(c(X_tem), probs = c(seq(0, 1, 1/(num_knots - 1)))) 
    }
    if(distri=="beta"){
      X_tem = array(rbeta(64*64*1000, 2,2),c(64,64,1000))
      knots = quantile(c(X_tem), probs = c(seq(0, 1, 1/(num_knots - 1)))) 
    }
    
    
    signal="BroadcasTR"
  }
  if(method=="TLR2"){
    signal="TLR-rescaled"
  }
  
  
  for(signal_i in 1:5){
    iditer= which(order(res_mat[,signal_i])==26)
    
    
    for(iter in iditer){
      
      plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
      
      if(method=="BNTR"){
        #BB=BB_l2norm_squre(result[[iter]][[signal_i]],knots,dimension = 2,original=1)-integration_of_entry_fun(result[[iter]][[signal_i]],knots)^2
        BB= fhatnorm_ten(result[[iter]][[signal_i]],knots)
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
#dev.off()