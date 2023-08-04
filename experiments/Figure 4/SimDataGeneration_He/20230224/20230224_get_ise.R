#library(BNTR)


source('./SimDataGeneration_He/20230224_get_inputs.R')
source('./SimDataGeneration_He/20230224_get_BB.R')
source('./SimDataGeneration_He/20230224_get_function.R')



# # 之前算ISE的方法
# BB = result[[iter]][[signal_i]]
# b0 = result[[iter]][[signal_i+5]]
# 
# #the entry is \Vert demean(m_j) - demean(\hat{m}_j) \Vert^2
# fun_diff_ten = functionpart_diff_ten(BB, case=signal_i, precision=1000, knots = knots, order = 4)
# #the entry is \int \hat{m}_j  
# fhat_intten = fhat_ten_int(BB, case=signal_i, precision = 1000,knots=knots, order = 4)
# #the entry is \int m_j (dose not include the intercept)
# ftrue_intten = ftrue_mean_ten(case=signal_i)
# # sum(fhat_intten)+b0 is \hat{intercept}, 1 - sum(ftrue_intten) is true intercept
# cons_diff = (sum(fhat_intten)+b0 - 1 - sum(ftrue_intten))^2
# # \sum_j \Vert demean(m_j) - demean(\hat{m}_j) \Vert^2
# fun_diff = sum(fun_diff_ten)
# # the entry is \int m_j^2 (dose not include the intercept) 
# ftrue_square_intten = ftrue_demean_square_int_ten(case=signal_i)
# # sum(ftrue_square_intten) is \sum \int m_j^2, sum(ftrue_intten) + 1 is the intercept
# ftrue_norm_squre = sum(ftrue_square_intten) + (sum(ftrue_intten) + 1)^2
# 
# true_b_squre = (sum(ftrue_intten) + 1)^2
# true_sum_mj = sum(ftrue_square_intten)
# ise_res = cons_diff + fun_diff
# ####




################### f hat part##################
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

fhat_function_ten <- function(X, BB,knots=c(0, 0.25, 0.5, 75, 1), order =4){
  tildeX = x_trans(X, knots = knots, order = order)
  res = array(0, dim(X))
  K = dim(tildeX)[length(dim(tildeX))]
  for(k in 1:K){
    res = res + tildeX[,,k] * BB[,,k]
  }
  return(res)
}

fhat_function_ten_linear <-function(X, BB){
  res = array(0, dim(X))
  res = X * BB
}

fhat_ten_int <- function(BB, precision = 10000,knots=c(0, 0.25, 0.5, 75, 1), order = 4){
  tem = 0
  p = dim(BB)[1:2]
  for(i in 1:precision/precision){
    X = array(i, p)
    tem = tem + fhat_function_ten(X, BB,knots=knots, order =order)
  }
  return(tem / precision)
}

fhat_ten_int_linear <- function(BB){
  return (0.5*BB)
}
################### f hat part##################


## BB[[i]] = BBprod[[i]][1,,]
ftrue_ten <- function(X, case, BBv){
  BBvv=list()
  for(i in 1:9){
    BBvv[[i]] = BBv[[i]][1,,]
  }
  BBv = BBvv
  
  if(case==1){
    res = BBv[[1]] * X
  }
  if(case==2){
    res = BBv[[2]] * f1(X) + BBv[[8]] * f2(X) + BBv[[9]] * f5(X)
  }
  if(case==3){
    res = BBv[[3]] * f2(X) + BBv[[4]] * f2(X) 
  }
  if(case==4){
    res = BBv[[5]] * f3(X) + BBv[[6]] * f5(X)
  }
  if(case==5){
    res = BBv[[7]] * f4(X)
  }
  
  
  # # for reference
  # y_all[[1]]=v[[1]]+t(ctprod(BBprod[[1]],X_data,2))
  # y_all[[2]]=v[[2]]+t(ctprod(BBprod[[2]],f1(X_data),2))+t(ctprod(BBprod[[8]],f2(X_data),2))+t(ctprod(BBprod[[9]],f5(X_data),2))
  # y_all[[3]]=v[[3]]+t(ctprod(BBprod[[3]],f2(X_data),2))+t(ctprod(BBprod[[4]],f2(X_data),2))
  # y_all[[4]]=v[[4]]+t(ctprod(BBprod[[5]],f3(X_data),2))+t(ctprod(BBprod[[6]],f5(X_data),2))
  # y_all[[5]]=v[[5]]+t(ctprod(BBprod[[7]],f4(X_data),2))
  # the dimension is from dim(X)
  
  return(res)
}

ftrue_ten_int <- function(case, BBv, precision = 100, p=c(32,32)){
  tem = 0 
  for(i in 1:precision/precision){
    X = array(i, p)
    tem = tem + ftrue_ten(X, case, BBv)
  }
  return(tem / precision)
}

ftrue_ten_int_squares <- function(case, BBv, precision = 100, p=c(32,32)){
  fture_mean = ftrue_ten_int(case, BBv, precision, p)
  
  tem = 0 
  for(i in 1:precision/precision){
    X = array(i, p)
    tem = tem + (ftrue_ten(X, case, BBv) - fture_mean)^2
  }
  return(tem / precision)
}

get_ftrue_ten_int_case_array <- function(BBv, case, p=c(32,32), precision=100){
  
  ftrue_ten_int_case_array = array(0,c(p, 5))
  #for(case in 1:5){
  ftrue_ten_int_case_array[,,case] = ftrue_ten_int(BBv=BBv, case=case, precision = precision, p=p)
  #}
  
  return(ftrue_ten_int_case_array)  
}




ftrue_demean_ten <- function(X, ftrue_ten_int_case_array, BBv, case=3){
  
  p = dim(X)[1:2]
  # make it mean zero
  res = ftrue_ten(X, case, BBv) - ftrue_ten_int_case_array[,,case]
  return(res)
}

#check if there exits problems
#p_test = c(64,64)
#ftrue_ten_int_case_array = get_ftrue_ten_int_case_array(p=p_test, precision=100)
#ftrue_demean_ten(X=matrix(rnorm(32*32),p_test), ftrue_ten_int_case_array, case=4)




final_ise <- function(b0, BB, v, BBv, case=1, precision=1000, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4){
  
  p = dim(BB)[1:2]
  fhat_intten = fhat_ten_int(BB,  precision = precision,knots=knots, order = order)
  ftrue_ten_int_case_array = get_ftrue_ten_int_case_array(case=case, BBv=BBv, p=p, precision=precision)
  
  # functionpaart_diff 
  tem = 0
  for (i in 1:precision/precision){
    X=array(i, p)
    Xtem = fhat_function_ten(X, BB,knots=knots, order =order) - fhat_intten
    #norm_ten = norm_ten + Xtem^2
    f_fun_hat = Xtem
    #f_fun_true =  ftrue_demean_ten(X, case)
    f_fun_true = ftrue_demean_ten(X, ftrue_ten_int_case_array, BBv, case=case)
    tem = tem +  (f_fun_hat -  f_fun_true)^2
  }
  #norm_ten = norm_ten/precision
  tem = tem/precision
  tem = tem
  #res = list(tem, norm_ten)
  
  fun_diff_ten = tem
  
  
  fun_diff = sum(fun_diff_ten)
  cons_diff = (sum(fhat_intten)+b0 - v - sum(ftrue_ten_int_case_array[,,case]))^2
  
  res = list(ise = fun_diff+cons_diff, fun_diff=fun_diff, cons_diff=cons_diff)
  
  
  return(res)
} 



final_ise_linear <- function(b0, BB, v, BBv, case=1, precision=1000, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4){
  
  p = dim(BB)[1:2]
  #fhat_intten = fhat_ten_int(BB,  precision = precision,knots=knots, order = order)
  fhat_intten = fhat_ten_int_linear(BB)
  ftrue_ten_int_case_array = get_ftrue_ten_int_case_array(case=case, BBv=BBv, p=p, precision=precision)
  
  # functionpaart_diff 
  tem = 0
  for (i in 1:precision/precision){
    X=array(i, p)
    #Xtem = fhat_function_ten(X, BB,knots=knots, order =order) - fhat_intten
    Xtem = fhat_function_ten_linear(X, BB) - fhat_intten
    #norm_ten = norm_ten + Xtem^2
    f_fun_hat = Xtem
    #f_fun_true =  ftrue_demean_ten(X, case)
    f_fun_true = ftrue_demean_ten(X, ftrue_ten_int_case_array, BBv, case=case)
    tem = tem +  (f_fun_hat -  f_fun_true)^2
  }
  #norm_ten = norm_ten/precision
  tem = tem/precision
  tem = tem
  #res = list(tem, norm_ten)
  
  fun_diff_ten = tem
  
  
  fun_diff = sum(fun_diff_ten)
  cons_diff = (sum(fhat_intten)+b0 - v - sum(ftrue_ten_int_case_array[,,case]))^2
  
  res = list(ise = fun_diff+cons_diff, fun_diff=fun_diff, cons_diff=cons_diff)
  
  
  return(res)
} 

