library(MASS)




#zhao的文章是对每一个样本都估计一个均值和方差
#新讨论，那就用极大似然原理来估计（zhao的文章并没有说明用什么方法估计）
#再次讨论，用所有样本来估计总的协方差矩阵。均值就用每个样本值来代替
VU_est <- function(Xmat){
  
  pm = dim(Xmat)[1]
  p_m = dim(Xmat)[2]
  if(length(dim(Xmat)) > 1){
    n = dim(Xmat)[3]
  }else{
    n=1
  }
  
  X_mean = matrix(0,pm,p_m)
  
  for(i in 1:pm){
    for(j in 1:p_m){
      X_mean[i,j] = mean(Xmat[i,j,])
    }
  }
  for(i in 1:n){
    Xmat[,,i] = Xmat[,,i] - X_mean
  }
  
  set.seed(2008)
  V0 = matrix(rnorm(p_m*p_m),p_m,p_m)
  U = 0
  U_new = 0
  
  V=V0
  
  for(iter in 1:100){
    Vinv = ginv(V)
    for(i in 1:n){
      U_new = U_new + Xmat[,,i]%*%Vinv%*%t(Xmat[,,i])
    }
    U_new = U_new/(p_m*n)
    Uinv = ginv(U_new)
    V_new = 0
    for(i in 1:n){
      V_new = V_new + t(Xmat[,,i])%*%Uinv%*%Xmat[,,i]
    }
    V_new = V_new/(pm*n)
    
    epsilonU = sum((U_new - U)^2)
    epsilonV = sum((V_new - V)^2)
    print(c(epsilonU, epsilonV))
    if(epsilonU < 0.01 && epsilonV < 0.01){
      V = V_new
      U = U_new
      #print("yes")
      break
    }
    V = V_new
    U = U_new
  }
  
  m_vec = c(X_mean)
  Cov = kronecker(V,U)
  res = list(m_vec = m_vec, Cov = Cov)
  return(res)
}
#

#symmetric Kullback-Leibler (KL) diver- gence
#Sigma_inv = ginv(Sigma^{-1})
sym_KL <- function(mu,Sigma_inv, Xmat1_vec, Xmat2_vec){
  #term1 = Xmat1_vec - mu
  #term2 = Xmat2_vec - mu
  #res = 0.5*term1%*%Sigma_inv%*%term1 +0.5*term2%*%Sigma_inv%*%term2 - 2*0.5 * length(mu)
  #try gaussian kernal
  #res  = sum((Xmat1_vec- Xmat2_vec)^2)
  #因为我们用的统一Sigma_inv,均值的话用自身来当估计
  term = Xmat1_vec - Xmat2_vec
  res =  t(term)%*%Sigma_inv%*%term
  
  return(res)
}


n_train = 10
n_pre = 2
p_use = c(2,2)
set.seed(100)
y = rnorm(n_train)
Xmat = array(runif(prod(c(p_use,n_train))),c(p_use,n_train))
Xmat_list = list()
Xmat_list[[1]] = Xmat
Xmat_list[[2]] = Xmat

Xmat_pre = array(runif(prod(c(p_use,n_pre))),c(p_use,n_pre))
Xmat_list_pre = list()
Xmat_list_pre[[1]] = Xmat_pre
Xmat_list_pre[[2]] = Xmat_pre


library(rTensor)
for(i in 1:n_train){
  tem = as.tensor(Xmat[,,i])
  Xmat_list[[1]][,,i]=k_unfold(tem,1)@data
  Xmat_list[[2]][,,i]=k_unfold(tem,2)@data
  y[i] = 10*sum(Xmat_list[[1]][,,i])
}
y_pre=rnorm(n_pre)
for(i in 1:n_pre){
  tem = as.tensor(Xmat_pre[,,i])
  Xmat_list_pre[[1]][,,i]=k_unfold(tem,1)@data
  Xmat_list_pre[[2]][,,i]=k_unfold(tem,2)@data
  y_pre[i] = 10*sum(Xmat_list_pre[[1]][,,i])
}


#为了加速运算，不要一个一个的算KL距离，算整体的距离
KXX_base_fun <- function(Xmat_list, theta){
  
  D = length(theta)-1
  n = dim(Xmat_list[[1]])[3]
  KXX_base = array(0, c(n, n, D))
  res = list()
  uuu = list()
  cov_inverse = list()
  for(d in 1:D){
    res[[d]] = VU_est(Xmat_list[[d]]) 
    uuu[[d]] = res[[d]]$m_vec
    cov_inverse[[d]] = chol2inv(chol(res[[d]]$Cov)) #ginv(res[[d]]$Cov) # 这里有个求逆，很费时间
  }
  
  ppp = prod(dim(Xmat_list[[1]][,,1]))
  
  for(d in 1:D){
    temmatrix = matrix(0,(1+n)/2 *n, ppp)
    for(i in 1:n){
      for(j in 1:i){
        temmatrix[(i-1+1)/2*(i-1)+j,] = c(Xmat_list[[d]][,,i]) - c(Xmat_list[[d]][,,j])
        #print(c(i,j))
        #KXX_base[i,j,d] = 2*kldiv(Xmat[,,i],Xmat[,,j])
        #KXX_base[i,j,d] = sym_KL(uuu[[d]], cov_inverse[[d]],c(Xmat_list[[d]][,,i]),c(Xmat_list[[d]][,,j]))
      }
    }
      print("over")
      temmatrix2 = crossprod(t(temmatrix), cov_inverse[[d]]) * temmatrix
    for(i in 1:n){
        for(j in 1:i){
          print(c(i,j))
          tem = sum(temmatrix[(i-1+1)/2*(i-1)+j,])
          KXX_base[i, j, d] = tem
          KXX_base[j, i, d] = tem
        }
      }

  }
  return(KXX_base) 
}



KXX_base_fun_pre <- function(Xmat_list_pre, Xmat_list, theta){
  
  D = length(theta)-1
  n_pre = dim(Xmat_list_pre[[1]])[3]
  n = dim(Xmat_list[[1]])[3]
  KXX_base_pre = array(0, c(n_pre, n, D))
  res = list()
  uuu = list()
  cov_inverse = list()
  for(d in 1:D){
    res[[d]] = VU_est(Xmat_list[[d]]) 
    uuu[[d]] = res[[d]]$m_vec
    cov_inverse[[d]] = chol2inv(chol(res[[d]]$Cov)) #ginv(res[[d]]$Cov) # 这里有个求逆，很费时间
  }
  
  
  ppp = prod(dim(Xmat_list[[1]][,,1]))
  
  for(d in 1:D){
    temmatrix = matrix(0,n_pre *n, ppp)
    for(i in 1:n_pre){
      for(j in 1:n){
        temmatrix[(i-1+1)/2*(i-1)+j,] = c(Xmat_list_pre[[d]][,,i]) - c(Xmat_list[[d]][,,j])
        #print(c(i,j))
        #KXX_base[i,j,d] = 2*kldiv(Xmat[,,i],Xmat[,,j])
        #KXX_base[i,j,d] = sym_KL(uuu[[d]], cov_inverse[[d]],c(Xmat_list[[d]][,,i]),c(Xmat_list[[d]][,,j]))
      }
    }
    print("over")
    temmatrix2 = crossprod(t(temmatrix), cov_inverse[[d]]) * temmatrix
    for(i in 1:n_pre){
      for(j in 1:n){
        print(c(i,j))
        tem = sum(temmatrix[(i-1+1)/2*(i-1)+j,])
        KXX_base_pre[i, j, d] = tem
      }
    }
    
  }
  return(KXX_base_pre) 
}


kernal_fun <- function(Xmat_list,theta, KXX_base){
  alpha = theta[1]
  betam = theta[-1]
  D = length(theta)-1
  
  #exp 里面的分子 ，包括-1/2
  #KXX_base = KXX_base_fun(Xmat_list,theta)
  n = dim(Xmat_list[[1]])[3]
  KXX_final = matrix(0,n,n)
  
  for(d in 1:D){
    KXX_final =   KXX_final  -0.5*KXX_base[,,d]/betam[d]^2
  }
  KXX_final = alpha^(2*D)*exp(KXX_final)
  return(KXX_final)
}


kernal_fun_pre <- function(Xmat_list_pre,Xmat_list,theta, KXX_base_pre){
  alpha = theta[1]
  betam = theta[-1]
  D = length(theta)-1
  
  #exp 里面的分子 ，包括-1/2
  #KXX_base = KXX_base_fun(Xmat_list,theta)
  n_pre = dim(Xmat_list_pre[[1]])[3]
  n = dim(Xmat_list[[1]])[3]
  KXX_final = matrix(0,n_pre,n)
  
  for(d in 1:D){
    KXX_final =   KXX_final  -0.5*KXX_base_pre[,,d]/betam[d]^2
  }
  KXX_final = alpha^(2*D)*exp(KXX_final)
  return(KXX_final)
}


#theta的导数，不包含sigma 
grad_kernal_fun<-function(Xmat_list,theta,KXX_base){
  alpha = theta[1]
  betam = theta[-1]
  D = length(theta)-1
  
  #alpha_all =alpha^(2*D)
  #exp 里面的分子, 不包括-1/2betam
  #D+1的原因是, 1个alpha, D个betam
  #alpha 在第一个位置，betam在后面的D个位置
  #KXX_base = KXX_base_fun(Xmat_list,theta)
  n = dim(Xmat_list[[1]])[3]
  KXX_grad = array(0,c(n,n,D+1))
  KXX_final = matrix(0,n,n)
  for(d in 1:D){
    KXX_final =   KXX_final - 0.5*KXX_base[,,d]/betam[d]^2
  }
  #KXX_final = \sum_{m=1}^M (-\frac{Cm}{2\bbeta_m^2})
  
  KXX_grad[,,1] = 2*D * alpha^(2*D-1) * exp(KXX_final)
  
  for(d in 1:D){
    KXX_grad[,,d+1] = alpha^(2*D)*exp(KXX_final)*(KXX_base[,,d]/betam[d]^3 ) 
  }
  
  #KXX_grad[,,D+2]=diag(n)
  return(KXX_grad)
}

grad_par <- function(Xmat_list,Xmat_list_pre,theta,y,y_pre,sigma,KXX_base,KXX_base_pre){
  
  n = length(y)
  grad_mat = grad_kernal_fun(Xmat_list,theta,KXX_base)
  Ky = kernal_fun(Xmat_list,theta,KXX_base) 
  Ky_inv = ginv( Ky + sigma^2 * diag(n))
  
  
  Kysigma = Ky + sigma^2*diag(n) 
  tem1= -2 * sigma*sum(diag(Ky_inv))+  2*sigma*sum(crossprod(y, Ky_inv)^2)
  #tem2=-sum(diag(2 * sigma *diag(n) %*% ginv(Kysigma)))+t(y)%*%ginv(Kysigma)%*%(2*sigma*diag(n))%*%ginv(Kysigma)%*%y
  
  grad_theta = theta
  D = length(theta)-1
  for(d in 1:(D+1)){
    tem = crossprod(t(Ky_inv), grad_mat[,,d])
    #grad_theta[d] = t(y) %*% tem %*% Ky_inv %*% y - sum(diag(tem))
    grad_theta[d] =  crossprod(t(crossprod(y, tem)), Ky_inv) %*%y- sum(diag(tem))
      #t(y) %*% tem %*% Ky_inv %*% y - sum(diag(tem))
  }
  grad_theta = 0.5 * grad_theta
  Ky_pre = kernal_fun_pre(Xmat_list_pre,Xmat_list,theta,KXX_base_pre) 
  y_pre_k = crossprod(t(Ky_pre), Ky_inv) %*%y 
  print(y_pre_k)
  mspe = sum((y_pre - y_pre_k)^2)/length(y_pre)
  
  #print(c(grad_theta,tem1,obj_fun(y,sigma,Ky),mspe))
  return(c(grad_theta,tem1,obj_fun(y,sigma,Ky),mspe))
}



obj_fun <- function(y,sigma,Ky){
  n = length(y)
  Kysigma = Ky + sigma^2*diag(n)
  tem=determinant(Kysigma)
  obj = - 0.5*tem$modulus[1] - 0.5*t(y) %*% ginv(Kysigma) %*%y
  return(obj)
}


prediction_y <- function(y, y_pre, kXnew_pre, ky_sigma){
  
  
  y = KXnewX %*% ginv(KXX + sigman2) %*%y
  return(y)
}






#关于
# grad_descent <-function(){
#   learning_rate = 0.01
#   set.seed(1993)
#   theta_sigma= rep(1,4)
#   D2 = length(theta_sigma)
#   D = D2 -2
#   grad_theta = rep(0,D2-1)
#   
#   for(iter in 1:4){
#     
#     theta = theta_sigma[-4]
#     sigma = theta_sigma[4]
#     Ky = kernal_fun(Xmat_list,theta)
#     KXX_grad=grad_kernal_fun(Xmat_list,theta)
#     
#     # #让sigma 的梯度最为接近0, 似乎不对，这样的化，直接令sigma趋于0，其梯度就为0了，因此不这样做。
#     # grad_sigma = c()
#     # for(iii in 1:1000){
#     #   sigma=iii/1000
#     #   Ky_inverse = ginv(Ky + sigma^2*diag(n))
#     #   tem = 2*sigma*Ky_inverse
#     #      grad_sigma[iii] = (0.5 * t(y) %*%tem%*%Ky_inverse %*%y - 0.5 * sum(diag(tem)))
#     # }
#     # index = which.min(abs(grad_sigma))
#     # sigma = index/1000
#     for(sigma in 1:10){
#       sigma=sigma/1000
#       obj_old = obj_fun(y,sigma,Ky)
#       print(obj_old)
#     }
#     
#     
#     
#     
#     for(d in 1:(D+1)){
#       Ky = kernal_fun(Xmat_list,theta)
#       Ky_inverse = ginv(Ky + sigma^2*diag(n))
#       tem = Ky_inverse%*%KXX_grad[,,d] 
#       grad_theta[d] =  0.5 * t(y) %*%tem%*%Ky_inverse %*%y - 0.5 * sum(diag(tem))
#     }
#     
#     
#     #似乎一起更新的化，跑不动。所以先把sigma 给更新了
#     obj_old = obj_fun(y,sigma,Ky)
#     print(obj_old)
#     obj_new = c()
#     for(inter_iter  in 1:10){
#       theta_tem = theta - 2*0.5^inter_iter * grad_theta
#       Ky_tem = kernal_fun(Xmat_list,theta_tem)
#       obj_new[inter_iter] = obj_fun(y,sigma,Ky_tem)
#       if(obj_new >= obj_old + 0.25 * 2*0.5^inter_iter * sum(theta^2 )){
#         break
#       }
#       index = which.max(obj_new)
#       theta = theta - 2*0.5^index * grad_theta
#       
#     }
#     
#     #determinant(kernal_fun(Xmat_list,theta))
#     
#   }
#   
#   
#   
#   
# }





#Ky = Kyy

#数值结果表明，并不是凸函数
# for(i in 1:1000){
#   sigma = i*0.001
#   #sigma=0
#   n = length(y)
#   Kysigma = Ky + sigma^2*diag(n) 
#   #print(sum(diag(2 * sigma *diag(n) * ginv(Kysigma))))
#   tem1=-sum(diag(2 * sigma *diag(n) %*% ginv(Kysigma)))+t(y)%*%ginv(Kysigma)%*%(2*sigma*diag(n))%*%ginv(Kysigma)%*%y
#   tem2 = obj_fun(y,sigma,Ky) 
#   print(c(tem1,tem2))
# }


# grad_test <- function(){
#   
#   Ky = Kyy + sigma
#   n = length(y)
#   t(y)%*%ginv(Kyy+sigma)%*%(2*sigma*diag(n))%*%ginv(Kyy+sigma)%*%y
# }

# KXnewX_base_fun <- function(Xnew, Xmat, theta){
#   
#   D = length(theta)-1
#   n = dim(Xmat)[3]
#   n_new = dim(Xnew)[3]
#   KXX_base = array(0, c(n_new, n , D))
#   for(i in 1:n_new){
#     for(j in 1:n){
#       for(d in 1:D){
#         KXX_base[i,j,d] = 2*kldiv(Xnew[,,i],Xmat[,,j])
#       }
#     }
#   }
#   KnewX = KXX_function(theta,KXX_base)
#   return(KnewX) 
# }


# KXX_function = function(theta, KXX_base){
#   D = length(theta)-1
#  # n = dim(Xmat)[3]
#   alpha = theta[1]^(2*D)
#   beta = 2*theta[-1]^2
#   tem = 0
#   for(d in 1:D){
#   tem = tem - KXX_base[,,d]/beta[d]
#   }
#   tem = alpha *exp(tem)
#   return(tem)
# }



# best_pre_y <- function(sigman2, y, Xnew, Xmat){
#   
# }

# 
# log_like = function(sigma,y){
#   sigmainv = ginv(sigma)
#   logterm = determinant(sigma,logarithm = TRUE)
#   n = length(y)
#   tem = 0
#   for(i in 1:n){
#     tem = tem + t(y) %*% sigmainv %*% y
#   }
#   res = - n * logterm - tem
#   return(res)
# }
# 
# hyper_par = function(a){
#   return(0)
# }




