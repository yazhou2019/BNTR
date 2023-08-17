load("~/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/B1.Rdata")
load("~/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/B2.Rdata")
load("~/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/B3.Rdata")
load("~/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/B41.Rdata")
load("~/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/B42.Rdata")
load("~/Desktop/JRSSBreview2/upload_JRSSB_2nd_respond/butterfly.Rdata")

CASE1 = B1
CASE2 = B2
CASE3 = B3
CASE4 = B41+B42
CASE5 = butterfly

fhatnorm_ten <- function(BB, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4) {
  pk = dim(BB)
  fhatsqure_intten <- apply(BB, c(1:(length(pk) - 1)), fhatsqure_int, knots = knots, order = order)
  fhat_intten <- apply(BB, c(1:(length(pk) - 1)), fhat_int, knots = knots, order = order)
  norm_ten <- sqrt(fhatsqure_intten - 2 * fhat_intten^2 +  fhat_intten^2)
  return(norm_ten)
}


x_trans <- function(X_sample, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4) {
  
  num_knots <- length(knots)
  
  # number of basis
  K <- num_knots - 2 + order - 1
  if (K == 1)
    return(X_sample)
  
  p <- dim(X_sample)
  BX <- array(0, c(p, K))
  str_slice <- c()
  for (j in 1:(length(p)+1)) {
    if (j == (length(p)+1)) {
      str_slice <- c(str_slice, "k")
    } else {
      str_slice <- c(str_slice, ",")
    }
  }
  str_slice <- paste(str_slice, collapse = " ")
  str_text1 <- paste("BX[", str_slice, "]=X_sample^k", collapse = " ")
  str_text2 <- paste("BX[", str_slice, "]=((X_sample - knots[k-order+1])*((X_sample - knots[k-order+1])>0))^(order - 1)",
                     collapse = " ")
  
  for (k in 1:K) {
    if (k < order) {
      eval(parse(text = str_text1))
    } else {
      eval(parse(text = str_text2))
    }
  }
  return(BX)
}




fhat_function <- function(x, coef=rep(1,6), knots=c(0, 0.25, 0.5, 75, 1), order =4){
  tildex = x_trans(x, knots = knots, order = order)
  res = sum(coef*tildex)
  return(res)
}

fhat_function_ten <- function(X, BB,knots=c(0, 0.25, 0.5, 75, 1), order =4){
  tildeX = x_trans(X, knots = knots, order = order)
  res = array(0, dim(X))
  K = dim(tildeX)[length(dim(tildeX))]
  for(k in 1:K){
    res = res + tildeX[,,k] * BB[,,k]
  }
  return(res)
}

fhat_ten_int <- function(BB, case=2, precision = 10000,knots=c(0, 0.25, 0.5, 75, 1), order = 4){
  tem = 0
  for(i in 1:precision/precision){
    X = array(i, c(64,64))
    tem = tem + fhat_function_ten(X, BB,knots=knots, order =order)
  }
  return(tem / precision)
}

  


fhat_function_ten_n <- function(X, BB,knots=c(0, 0.25, 0.5, 75, 1), order =4){
  tildeX_n = tildePhiX_trans(X, knots = knots, order = order)
  res = array(0, dim(X))
  K = dim(BB)[length(dim(BB))]
  for(k in 1:K){
    res = res + tildeX[,,k,] * BB[,,k]
  }
  return(res)
}

start_time = Sys.time()
#tildeX_n = tildePhiX_trans(X, knots = knots, order = order)
#for(i in 1:10000){fhat_function_ten(X, BB)}
#for(i in 1:10000000){fconst_demean(1, 0.5)}
end_time = Sys.time()
time_taken = end_time - start_time
print(time_taken)

fconst <- function(x){
  return(x)
}
f1 <- function(x){
  return(x + 0.6*sin(2*pi*(x-0.5)^2))
}
f2 <- function(x){
  return(x + 0.3*cos(2*pi*x))
}

f_int <- function(f, precision = 10000){
  s = 0
  for(i in 1:precision/precision){
    s = s + f(i)
  }
  return(s/precision)
}

fconst_int = 0.5
f1_int = f_int(f1,100000)
f2_int = f_int(f2,100000)

fconst_demean <- function(x){
  return(fconst(x) - fconst_int)
}
f1_demean <- function(x){
  return(f1(x) - f1_int)
}
f2_demean <- function(x){
  return(f2(x) - f2_int)
}

fconts_demean_square <- function(x){
  return(fconst_demean(x)^2)
}
f1_demean_square <- function(x){
  return(f1_demean(x)^2)
}
f2_demean_square <- function(x){
  return(f2_demean(x)^2)
}

fconts_demean_square_int = f_int(fconts_demean_square, 100000)
f1_demean_square_int =  f_int(f1_demean_square,100000)
f2_demean_square_int =  f_int(f2_demean_square,100000)


#B1, B2, B3, B41, B41, butterfly have been loaded in the environment
ftrue_demean <- function(x, i=3,j=4, case=3){
  # make it mean zero
  if(case==1){
    return(B1[i,j]* (fconst(x) - 0.5) )
  }
  if(case==2){
    return(B2[i,j]*(f1(x) - f1_int))
  }
  if(case==3){
    return(B3[i,j]*(f1(x) - f1_int))
  }
  if(case==4){
    return(B41[i,j]*(f1(x)-f1_int)+B42[i,j]*(f2(x)- f2_int) )
  }
  if(case==5){
    return(butterfly[i,j] *(f1(x) - f1_int))
  }
}

ftrue_demean_ten <- function(X, case=3){
  
  # make it mean zero
  if(case==1){
    return(B1 * (fconst(X) - 0.5) )
  }
  if(case==2){
    return(B2 * (f1(X) - f1_int))
  }
  if(case==3){
    return(B3 * (f1(X) - f1_int))
  }
  if(case==4){
    return(B41 * (f1(X)-f1_int)+B42 * (f2(X)- f2_int) )
  }
  if(case==5){
    return(butterfly * (f1(X) - f1_int))
  }
}

ftrue_mean_ten <- function(case=3){
  # make it mean zero
  if(case==1){
    return(B1 * (0.5) )
  }
  if(case==2){
    return(B2 * (f1_int))
  }
  if(case==3){
    return(B3 * (f1_int))
  }
  if(case==4){
    return(B41 * (f1_int)+B42 * (f2_int) )
  }
  if(case==5){
    return(butterfly * (f1_int))
  }
}




ftrue_demean_square_int_ten <- function(case =3){
  # make it mean zero
  if(case==1){
    return(B1 * fconts_demean_square_int )
  }
  if(case==2){
    return(B2 * f1_demean_square_int)
  }
  if(case==3){
    return(B3 * f1_demean_square_int)
  }
  if(case==4){
    return(B41 * f1_demean_square_int+B42 * f2_demean_square_int )
  }
  if(case==5){
    return(butterfly * f1_demean_square_int)
  }
}



# \int f_0(\hat f - \int \hat f) 

# inter_term_ten <- function(BB,case=1, precision=10000, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4){
#   pk = dim(BB)
#   tem1 = 0
#   tem2 = 0
#   tem3 = 0
#   #fhat_intten <- apply(BB, c(1:(length(pk) - 1)), fhat_int, knots = knots, order = order)
#   fhat_intten = fhat_ten_int(BB, case=case, precision = precision,knots=knots, order = order)
#   
#   for (i in 1:precision/precision){
#   X=array(i, c(64, 64))
#   Xtem = fhat_function_ten(X, BB,knots=knots, order =order) - fhat_intten
#   #tem2 = tem2 +ftrue_demean_ten(case) * Xtem 
#   tem3 = tem3 +  (Xtem - ftrue_demean_ten(X, case))^2
#   #tem1 = tem1 + Xtem^2
#   }
#   tem1 = tem1/precision
#   tem2 = tem2/precision
#   res = list(tem1, tem2)
#   return(res)
# }
  
abstract_ten <- function(BB, case=1, precision=10000, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4){

  fhat_intten = fhat_ten_int(BB, case=case, precision = precision,knots=knots, order = order)
  #fhat_intten <- apply(BB, c(1:(length(pk) - 1)), fhat_int, knots = knots, order = order)
  tem = 0
  #norm_ten = 0
  for (i in 1:precision/precision){
  X=array(i, c(64, 64))
  Xtem = fhat_function_ten(X, BB,knots=knots, order =order) - fhat_intten
  #norm_ten = norm_ten + Xtem^2
  tem = tem +  (Xtem- ftrue_demean_ten(X, case))^2
  
  }
  #norm_ten = norm_ten/precision
  tem = tem/precision
  tem = tem^0.5 
  #res = list(tem, norm_ten)
  return(tem)
} 



#mat = array(1:64^2,c(64,64))
ordermap2tensor<-function(x){
  j = as.integer(x/64)
  i = x - 64*j
  return(c(i, j+1))
}
  
  

library(ggplot2)

# coef = BB[i, j]
ftrue_and_hat_plot <-function(BB, i, j, case=2, knots=c(0, 0.25, 0.5, 0.75, 1)){
  
  y_hat <- c()
  y_true <- c()
  precision = 100
  
  hat_minus_true = c()
  
  for(index in 1:precision){
  x = index/precision
  y_hat[index] = fhat_function(x, BB[i,j,],knots=knots) - fhat_int(BB[i,j,], knots= knots)
  y_true[index] = ftrue_demean(x, i=i,j=j, case=case)
  hat_minus_true[index] = abs(y_hat[index] -  y_true[index])
  }
  
  print((sum(hat_minus_true^2)*1/precision)^0.5)
  
  iteration = c(1:precision/precision,1:precision/precision)
  ob_value = c(y_true, y_hat)
  
  group = c(rep('true', precision), rep('esti', precision))
  
  df = data.frame(iteration , ob_value, group)
  
  res = ggplot(df,aes(x=iteration, y=ob_value,colour=group, group=group))+geom_line()+labs(x="x",y="function value") + scale_color_manual(values=c("red","blue"))
  res = res + ylim(-1,1)
  # plot(1:precision, y_true)
  # lines(1:precision, y_hat, type = "s")
  # print("summation of y_true")
  # print(sum(y_true))
  # print("summation of y_hat")
  # print(sum(y_hat))
  return(res)
}

  
  
histogram_densitycheck <-function(X){
  X =data.frame(X=c(X))
  p <- ggplot(X, aes(x=X))  + geom_histogram(breaks=seq(0,1,1/20)) + xlim(0,1)
  p
}


library(stringr)

plot_fhatminusftrue <-function(signal_i=1, iter=1, qseq = c(0.95, 0.975, 0.99, 1), result=NA, knots = NA, n_data=1000, onlyone=0,BBB=NA){
  p=c(64,64)
  
  if(is.na(sum(knots))){
  if(n_data==1000){
  num_knots=6
  }
  if(n_data==750||n_data==500){
  num_knots=5
  }
  if(n_data==250){
  num_knots=4  
  }  
    
    
    if(distri == "tnorm"){
      X = array(rtnorm(64*64*1000, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,1000))
      knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
    }
    if(distri == "unif"){
      X = array(runif(64*64*1000), c(64,64,1000))
      knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
    }
    
  #X = array(rtnorm(64*64*1000, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,1000))
  #knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
  }
  
  #signal_i = 5
  plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
  
  ## the previos method to calculate norm matrix
  #BB <- fhatnorm_ten(result[[2]][[signal_i]], knots=c(0,0.2,0.4,0.6,0.8,1))
  #rasterImage(1-BB, 1, 1, p[1], p[2]) 
  
  
  if(is.na(sum(BBB))==1){
  abstract_tensornorm = abstract_ten(result[[iter]][[signal_i]], precision=1000, case=signal_i, knots = knots)
  }else{
  abstract_tensornorm = abstract_ten(BBB, precision=1000, case=signal_i, knots = knots)
  }
  #rasterImage(1-abstract_tensornorm, 1, 1, p[1], p[2]) 
  #rasterImage(abstract_tensornorm, 1, 1, p[1], p[2]) 
  print(str_c("sum all error=",sum(abstract_tensornorm)))
  
  if(onlyone==0){
  use_index_se = order(abstract_tensornorm)[as.integer(stats::quantile(1:64^2, probs=qseq))]
  }
  if(onlyone==1){
  if(signal_i==1){caseuse = CASE1}
  if(signal_i==2){caseuse = CASE2}
  if(signal_i==3){caseuse = CASE3}
  if(signal_i==4){caseuse = CASE4}
  if(signal_i==5){caseuse = CASE5}
  use_index_se = order(abstract_tensornorm[caseuse==1])[as.integer(stats::quantile(1:length(caseuse[caseuse==1]), probs=qseq))]  
  use_index_se = which(caseuse==1)[use_index_se]
  print(stats::quantile(c(abstract_tensornorm[caseuse==1]), probs=qseq))
  }
  
  res = list()
  for(i in 1:4){
    #abstract_tensornorm[use_index_se[5]]
    index = ordermap2tensor(use_index_se[i])
    print("")
    print(c(abstract_tensornorm[index[1],index[2]], index[1],index[2]))
    #res[[i]] = ftrue_and_hat_plot(result[[iter]][[signal_i]], index[1], index[2], case=signal_i, knots=knots)
    if(is.na(sum(BBB))==1){
      res[[i]] = ftrue_and_hat_plot(result[[iter]][[signal_i]], index[1], index[2], case=signal_i, knots=knots)
    }else{
      res[[i]] = ftrue_and_hat_plot(BBB, index[1], index[2], case=signal_i, knots=knots)
    }
    
    #titlename =  str_c("value=",  round(abstract_tensornorm[index[1],index[2]], 3))
    titlename = str_c("quantilte=", qseq[i])
    res[[i]] = res[[i]] + labs(title = titlename)
  } 
  library(Rmisc)
  #return(multiplot(res[[1]],res[[2]],res[[3]],res[[4]],cols=2))
  return(res)
}

plot_fhatminusftrue_vary_sample <-function(signal_i=1, iter=1, 
                                           qseq = c(0.95, 0.975, 0.99, 1), 
                                           result=NA, knots = NA, n_data=1000, 
                                           onlyone=0,BBB=NA,distri="unif",
                                           index=NA
                                           ){
  p=c(64,64)
  
  if(is.na(sum(knots))){
    if(n_data==1000){
      num_knots=6
    }
    if(n_data==750||n_data==500){
      num_knots=5
    }
    if(n_data==250){
      num_knots=4  
    }  
    
    
    if(distri == "tnorm"){
      X = array(rtnorm(64*64*1000, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,1000))
      knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
    }
    if(distri == "unif"){
      X = array(runif(64*64*1000), c(64,64,1000))
      knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
    }
    
    #X = array(rtnorm(64*64*1000, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,1000))
    #knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
  }
  
  #plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
  
  if(is.na(sum(BBB))==1){
    abstract_tensornorm = abstract_ten(result[[iter]][[signal_i]], precision=1000, case=signal_i, knots = knots)
  }else{
    abstract_tensornorm = abstract_ten(BBB, precision=1000, case=signal_i, knots = knots)
  }
  print(str_c("sum all error=",sum(abstract_tensornorm)))
  
  if(onlyone==0){
    use_index_se = order(abstract_tensornorm)[as.integer(stats::quantile(1:64^2, probs=qseq))]
  }
  if(onlyone==1){
    if(signal_i==1){caseuse = CASE1}
    if(signal_i==2){caseuse = CASE2}
    if(signal_i==3){caseuse = CASE3}
    if(signal_i==4){caseuse = CASE4}
    if(signal_i==5){caseuse = CASE5}
    use_index_se = order(abstract_tensornorm[caseuse==1])[as.integer(stats::quantile(1:length(caseuse[caseuse==1]), probs=qseq))]  
    use_index_se = which(caseuse==1)[use_index_se]
    print(stats::quantile(c(abstract_tensornorm[caseuse==1]), probs=qseq))
  }
  
  res = list()
  for(i in 1:length(qseq)){
    #abstract_tensornorm[use_index_se[5]]
    if(is.na(index)==1){
    index = ordermap2tensor(use_index_se[i])
    }
    print("")
    print(c(abstract_tensornorm[index[1],index[2]], index[1],index[2]))
    if(is.na(sum(BBB))==1){
      res[[i]] = ftrue_and_hat_plot(result[[iter]][[signal_i]], index[1], index[2], case=signal_i, knots=knots)
    }else{
      res[[i]] = ftrue_and_hat_plot(BBB, index[1], index[2], case=signal_i, knots=knots)
    }
    
    #titlename =  str_c("value=",  round(abstract_tensornorm[index[1],index[2]], 3))
    titlename = str_c("quantilte=", qseq[i])
    res[[i]] = res[[i]] + labs(title = titlename)
  } 
  #library(Rmisc)
  #return(multiplot(res[[1]],res[[2]],res[[3]],res[[4]],cols=2))
  final_res = list(res=res, index=index)
  #return(res)
  return(final_res)
}

x = matrix(rep(1:64, times=64), 64,64)
y = matrix(rep(1:64, each=64), 64,64)

#res=list()
#res[[1]]=poly.image(x,y,BB, zlim = c(0,1),col = redrangehex(200,3), xlab=NA, ylab=NA,xaxt = "n", yaxt="n")
#res[[2]] = res[[1]]


#data = data.frame(x = c(x), y = c(y), z=c(BB))

#ggplot() + geom_raster(data = data , aes(x = x, y = y, fill = z))  + scale_color_manual(values= redrangehex(200,3))
# only red

redrangehex <- function(precision, coloruse = 1){
rgp2hexcode <- function(x, coloruse = 1){
  
x_use = rep(0,3)
x_use[coloruse] = 255
x_use[-coloruse] = x
  
return(rgb(x_use[1], x_use[2], x_use[3], maxColorValue=255))
}
res = c()
colrange = seq(0,255,255/precision)
colrange = colrange[length(colrange):1]

for(i in 1:length(colrange)){
  res[i] = rgp2hexcode(colrange[i], coloruse)
}
return(res)
}


# try case 3
m_true <-function(xvec){
  X_data = array(xvec, c(64, 64,1))
  X_data_reg = X_data+adjust_nonlinear*0.5*cos(2*pi*X_data)
return(1+t(ctprod(BBprod[[3]], X_data_reg,2))+t(ctprod(BBprod[[4]],X_data_reg,2)))
}

m_test <-function(xvec){
  return(sum(xvec))
}



functionpart_diff_ten <- function(BB, case=1, precision=1000, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4){
  
  fhat_intten = fhat_ten_int(BB, case=case, precision = precision,knots=knots, order = order)
  #fhat_intten <- apply(BB, c(1:(length(pk) - 1)), fhat_int, knots = knots, order = order)
  tem = 0
  for (i in 1:precision/precision){
    X=array(i, c(64, 64))
    Xtem = fhat_function_ten(X, BB,knots=knots, order =order) - fhat_intten
    #norm_ten = norm_ten + Xtem^2
    f_fun_hat = Xtem
    f_fun_true =  ftrue_demean_ten(X, case)
    tem = tem +  (f_fun_hat -  f_fun_true)^2
    
  }
  #norm_ten = norm_ten/precision
  tem = tem/precision
  tem = tem
  #res = list(tem, norm_ten)
  return(tem)
} 


two_est_diff_ten <- function(BB1, BB2, b1,b2, precision=1000, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4){
  
  fhat_intten1 = fhat_ten_int(BB1, precision = precision,knots=knots, order = order)
  fhat_intten2 = fhat_ten_int(BB2, precision = precision,knots=knots, order = order)
 
  tem = 0
  for (i in 1:precision/precision){
    X=array(i, c(64, 64))
    f_fun_hat1 = fhat_function_ten(X, BB1,knots=knots, order =order) - fhat_intten1
    f_fun_hat2 = fhat_function_ten(X, BB2,knots=knots, order =order) - fhat_intten2
    tem = tem +  (f_fun_hat1 -  f_fun_hat2)^2
    
  }
  tem = tem/precision
  functionpart = sum(tem)
  constantpart = (sum(fhat_intten1) + b1 - sum(fhat_intten2) - b2)^2
  ISE = functionpart + constantpart
  res = list(ISE=ISE, functionpart = functionpart, constantpart=constantpart)
  return(res)
} 



plot_fhatminusftrue_all <-function(signal_i=1, iter=1, result=NA, knots = NA, n_data=1000, onlyone=1,BBB=NA, distri="tnorm"){
  p=c(64,64)
  
  if(is.na(sum(knots))){
    if(n_data==1000){
      num_knots=6
    }
    if(n_data==750||n_data==500){
      num_knots=5
    }
    if(n_data==250){
      num_knots=4  
    }  
    
    if(distri == "tnorm"){
    X = array(rtnorm(64*64*1000, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,1000))
    knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
    }
    if(distri == "unif"){
    X = array(runif(64*64*1000), c(64,64,1000))
    knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
    }
    
  }
  
  #signal_i = 5
  plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")
  
  
  if(is.na(sum(BBB))==1){
    abstract_tensornorm = abstract_ten(result[[iter]][[signal_i]], precision=1000, case=signal_i, knots = knots)
  }else{
    abstract_tensornorm = abstract_ten(BBB, precision=1000, case=signal_i, knots = knots)
  }
  
  print(str_c("sum all error=",sum(abstract_tensornorm^2)))
  
  if(onlyone==0){
    use_index_se = order(abstract_tensornorm)[as.integer(stats::quantile(1:64^2, probs=qseq))]
  }
  if(onlyone==1){
    if(signal_i==1){caseuse = CASE1}
    if(signal_i==2){caseuse = CASE2}
    if(signal_i==3){caseuse = CASE3}
    if(signal_i==4){caseuse = CASE4}
    if(signal_i==5){caseuse = CASE5}
    use_index_se = order(abstract_tensornorm[caseuse==1]) #[as.integer(stats::quantile(1:length(caseuse[caseuse==1]), probs=qseq))]  
    use_index_se = which(caseuse==1)[use_index_se]
    #print(stats::quantile(c(abstract_tensornorm[caseuse==1]), probs=qseq))
  }
  
  res = list()
  res_index = list()
  for(i in 1:length(use_index_se)){
    index = ordermap2tensor(use_index_se[i])
    res_index[[i]] =index
    print("")
    print(c(abstract_tensornorm[index[1],index[2]], index[1],index[2]))
    if(is.na(sum(BBB))==1){
    res[[i]] = ftrue_and_hat_plot(result[[iter]][[signal_i]], index[1], index[2], case=signal_i, knots=knots)
    }else{
    res[[i]] = ftrue_and_hat_plot(BBB, index[1], index[2], case=signal_i, knots=knots)
    }
    #titlename =  str_c("value=",  round(abstract_tensornorm[index[1],index[2]], 3))
    tem_error = round(abstract_tensornorm[index[1],index[2]]^2,4)
    titlename = str_c('order=',i,',error=',tem_error, ",(",index[1],',',index[2],")" )
    res[[i]] = res[[i]] + labs(title = titlename)
  } 
  library(Rmisc)
  #multiplot(res[[1]],res[[2]],res[[3]],res[[4]],cols=2)
  res_final = list(plot =res, index=res_index)
  return(res_final)
}









abstract_ten_linear <- function(BB, case=1, precision=10000){
  tem = 0
  for (i in 1:precision/precision){
    X=array(i, c(64, 64))
    Xtem = BB*X - BB/2
    tem = tem +  (Xtem- ftrue_demean_ten(X, case))^2
  }
  tem = tem/precision
  tem = tem^0.5 
  return(tem)
} 


ftrue_and_hat_plot_linear <-function(BB, i, j, case=2){
  
  y_hat <- c()
  y_true <- c()
  precision = 100
  
  hat_minus_true = c()
  
  for(index in 1:precision){
    x = index/precision
    y_hat[index] = x*BB[i,j] - BB[i,j]/2 
    y_true[index] = ftrue_demean(x, i=i,j=j, case=case)
    hat_minus_true[index] = abs(y_hat[index] -  y_true[index])
  }
  
  print((sum(hat_minus_true^2)*1/precision)^0.5)
  
  iteration = c(1:precision/precision,1:precision/precision)
  ob_value = c(y_true, y_hat)
  
  group = c(rep('true', precision), rep('esti', precision))
  
  df = data.frame(iteration , ob_value, group)
  
  res = ggplot(df,aes(x=iteration, y=ob_value,colour=group, group=group))+geom_line()+labs(x="x",y="function value") + scale_color_manual(values=c("red","blue"))
  res = res + ylim(-1,1)
  return(res)
}



plot_fhatminusftrue_vary_sample_linear <-function(
  sample_size_list,
  signal_i=1, 
  qseq = c(0.95, 0.975, 0.99, 1), 
  size_i=1, 
  method_i=1,
  onlyone=0,
  distri="unif",
  index = NA){
  p=c(64,64)
  
  
  BBB = sample_size_list[[size_i]][[method_i]][[signal_i]]
  
  abstract_tensornorm = abstract_ten_linear(BBB, precision=1000, case=signal_i)
  
  print(str_c("sum all error=",sum(abstract_tensornorm)))
  
  if(onlyone==0){
    use_index_se = order(abstract_tensornorm)[as.integer(stats::quantile(1:64^2, probs=qseq))]
  }
  if(onlyone==1){
    if(signal_i==1){caseuse = CASE1}
    if(signal_i==2){caseuse = CASE2}
    if(signal_i==3){caseuse = CASE3}
    if(signal_i==4){caseuse = CASE4}
    if(signal_i==5){caseuse = CASE5}
    use_index_se = order(abstract_tensornorm[caseuse==1])[as.integer(stats::quantile(1:length(caseuse[caseuse==1]), probs=qseq))]  
    use_index_se = which(caseuse==1)[use_index_se]
    print(stats::quantile(c(abstract_tensornorm[caseuse==1]), probs=qseq))
  }
  
  res = list()
  for(i in 1:length(qseq)){
    #abstract_tensornorm[use_index_se[5]]
    if(is.na(index)==1){
    index = ordermap2tensor(use_index_se[i])
    }
    print("")
    print(c(abstract_tensornorm[index[1],index[2]], index[1],index[2]))
    
    
    res[[i]] = ftrue_and_hat_plot_linear(BBB, index[1], index[2], case=signal_i)
    
    #titlename =  str_c("value=",  round(abstract_tensornorm[index[1],index[2]], 3))
    titlename = str_c("quantilte=", qseq[i])
    res[[i]] = res[[i]] + labs(title = titlename)
  } 
  res_final = list(res = res,index=index)
  return(res_final)
}



