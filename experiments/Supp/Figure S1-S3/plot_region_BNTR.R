#plot vary size



##plot for BNTR in various  sample size
p=c(64,64)
par(mar = c(0,0,0,0),mfcol=c(5,3),mai=c(0.05,0.25,0.3,0.2))


n_all=c(250,500,750,1000)
BB_mat <- array(0, c(64,64,4,3))
largest_value = matrix(0, 4, 4)
res_mat_f = matrix(0, 50, 5)

color_used = c(5/255,120/255,5/255)
scale_par = 0.8
for(n_i in 2:4){
  method="BNTR" # nonlinear model
  #method="TLR2" # linear model with rescaling
  n=n_all[n_i]
  n_size=str_c("n=",n)
  
  
  # if(n==250){
  #   if(method=='BNTR'){
  #     #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/SimResults/BNTR250_new_K4.Rdata")
  #     #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/MSE_ISE/res_mat250K4uniform.Rdata")
  #     
  #     load("./SimResults/BNTR250_new_K4.Rdata")
  #     load("./MSE_ISE/res_mat250K4uniform.Rdata")
  #     result_1_4 =result 
  #     res_mat_f[,1:4] = res_mat[,1:4]
  #     #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/MSE_ISE/res_mat250K4butterfly.Rdata")
  #     #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/SimResults/BNTR250_butterfly_K4.Rdata")
  #     load("./MSE_ISE/res_mat250K4butterfly.Rdata")
  #     load("./SimResults/BNTR250_butterfly_K4.Rdata")
  #     
  #     result_5 = result
  #     res_mat_f[,5] = res_mat[,5]}
  # }
  if(n ==500 || n==750){
    library(stringr)
    library(BNTR)
    source('./CollectISEs/ISE_comput_functions.R')
    
    path1=str_c("./SimResults/", method, "_", n, "_new.Rdata")
    path2=str_c("./SimResults/final_error_", method, "_", n, "_new.Rdata")
    
    load(path1)
    load(path2)
    result_1_4 =result 
    res_mat_f[,1:4] = final_error[1:4,]
    if(n==500){
      #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/SimResults/BNTR500_butterfly_K5.Rdata")
      #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/MSE_ISE/res_mat500K5butterfly.Rdata")
      load("./SimResults/BNTR500_butterfly_K5.Rdata")
      load("./MSE_ISE/res_mat500K5butterfly.Rdata")
      
    }
    if(n==750){
      #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/SimResults/BNTR750_butterfly_K5.Rdata")
      #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/MSE_ISE/res_mat750K5butterfly.Rdata")
      load("./SimResults/BNTR750_butterfly_K5.Rdata")
      load("./MSE_ISE/res_mat750K5butterfly.Rdata")
      
    }
    result_5 = result
    res_mat_f[,5] = res_mat[,5]
    
    
  }
  if(n==1000){
    load("./SimResults/BNTR1000_new_K7_2rd.Rdata")  
    load("./SimResults/res_mat1000K8.Rdata")
    
    result_1_4 =result 
    res_mat_f[,1:4] = res_mat[,1:4]
    
    #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/SimResults/BNTR1000_butterfly_K7.Rdata")
    #load("/Users/ya/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/NewResults/MSE_ISE/res_mat1000K7butterfly.Rdata")
    load("./SimResults/BNTR1000_butterfly_K7.Rdata")
    load("./MSE_ISE/res_mat1000K7butterfly.Rdata")
    result_5 = result
    res_mat_f[,5] = res_mat[,5]
    
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
    
    
    X_tem = array(runif(64*64*1000,0,1),c(64,64,1000))
    knots = quantile(c(X_tem), probs = c(seq(0, 1, 1/(num_knots - 1)))) 
    
    
    
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
        if(signal_i<=4){
          BB= fhatnorm_ten(result_1_4[[iter]][[signal_i]],knots)
          #BB_mat[,, signal_i, n_i]=BB
        }
        if(signal_i==5){
          BB= fhatnorm_ten(result_5[[iter]][[1]],knots)
        }
      } else {
        BB=result[[iter]][[signal_i]] 
        if(sum(BB)!=0){
          BB=sqrt(BB^2/12)
        }
      }
      
      #largest_value[signal_i, n_i]=max(BB)
      BBplot = array(1, c(64,64,3))
      BBplot[,,-2] = 1 - BB
      
      BBplot[,,1] = color_used[1]  * ((-1/color_used[1] + 1) *BB/scale_par + 1/color_used[1])
      BBplot[,,2] = color_used[2]  * ((-1/color_used[2] + 1) *BB/scale_par + 1/color_used[2]) 
      BBplot[,,3] = color_used[3] * ((-1/color_used[3] + 1) *BB/scale_par + 1/color_used[3])
      
      rasterImage(BBplot, 1, 1, p[1], p[2]) 
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
