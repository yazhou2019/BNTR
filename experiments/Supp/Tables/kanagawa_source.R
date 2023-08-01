library(inline)
library(microbenchmark)
Cprod = cxxfunction(signature(tm="NumericMatrix",
                              tm2="NumericMatrix"),
                    plugin="RcppEigen",
                    body="
                    NumericMatrix tm22(tm2);
                    NumericMatrix tmm(tm);
                    const Eigen::Map<Eigen::MatrixXd> ttm(as<Eigen::Map<Eigen::MatrixXd> >(tmm));
                    const Eigen::Map<Eigen::MatrixXd> ttm2(as<Eigen::Map<Eigen::MatrixXd> >(tm22));
                    Eigen::MatrixXd prod = ttm*ttm2;
                    return(wrap(prod)); ")


library(MASS)
library(matrixcalc)
#library(BNTR)
kanagawa_GP <- function(R=1, xkn, xkpre, y){
  N = 2000
  n = length(xkn)
  K = length(xkn[[1]])
  n_pre = length(xkpre)
  beta = 0.1 #在某个例子中，这个越小，拟合的越好
  f =  array(rnorm(R*K*n), c(R, K, n))
  fpre = array(0, c(R,K,n_pre))
  
  lambda = 1 #在某个例子中，这个越小，拟合的越好
  delta=10 
  #sampling
  # f[r,k,] =  mvrnorm(1,rep(0, 2),Sigma)
  Kmat_all = K_com_all(xkn,delta)/lambda
  Kmat_all_inv = Kmat_all
  for(k in 1:K){
    #Kmat_all_inv[,,k] = chol2inv(chol(Kmat_all[,,k]))
    Kmat_all_inv[,,k] = solve(Kmat_all[,,k])
  }
  Kmat_pre = K_com_pre(xkpre, xkn, delta)/lambda
  
  ypre = 0
  for(iter in 1:N){
    print(iter)
    for(k in 1:K){
      for(r in 1:R){
        #1, caculation a,b at first
        a = c()
        b = c()
        kvec = 1:K
        rvec = 1:R
        for(i in 1:n){
          #print(f[r,kvec[-k], i])
          a[i] = prod(f[r,kvec[-k], i])
          tem = 0
          for(r_tem in rvec[-r]){
            tem = tem +  prod(f[r_tem,,i]) 
          }
          b[i] = tem
        }
        # print(b)
        
        #since for each r, we use the same Krk
        Krk = Kmat_all[,,k]
        #print(a)
        #2，calculate the Cov and mean
        #by Woodbury matrix identity, we have the following result, Note that the formula in kanagawa's paper lacks an "inverse"
        Sigmark = chol2inv(chol(Kmat_all_inv[,,k] + diag(a^2)/beta))
        murk = Sigmark %*% (hadamard.prod(a,y-b))/beta
        #print(hadamard.prod(a,y-b))
        #print(murk)
        #3, sample frk
        f[r,k,] =  mvrnorm(1, murk, Sigmark)
      }
    }
    
    #prediction by posterier mean
    if(iter >=1000){
      tem = matrix(1, n_pre, R)
      for(k in 1:K){
        #f[,k,] R \times n
        #fre[,k,] R\times n
        fpre[,k,] = t(Cprod(Cprod(Kmat_pre[,,k], Kmat_all_inv[,,k]),t(matrix(c(f[,k,]),R,n))))
        #fpre[,k,] = t(Kmat_pre[,,k] %*% Kmat_all_inv[,,k] %*% t(matrix(c(f[,k,]),R,n)))
        tem = hadamard.prod(tem,t(fpre[,k,])) #n_pre \times R
      }
      ypre = ypre + rowSums(tem) 
    }
  }
  ypre = ypre/(N-1000)
  return(ypre)
}
#input 形式



kanagawa_GP_vali <- function(R_se, beta_se,lambda_se,delta_se, xkn_train, xkn_vali, xkn_test,y_train, y_vali, y_test){
  
  n_R_se = length(R_se)
  n_beta_se = length(beta_se)
  n_lambda_se = length(lambda_se)
  n_delta_se = length(delta_se) 
  n_tuning = n_R_se * n_beta_se * n_lambda_se * n_delta_se
  
  n_vali = length(y_vali)
  n_test = length(y_test)
  n_train = length(y_train)
  
  res_tem = list(y_vali = 0, y_test = 0)
  res_mat = matrix(0, n_tuning, 3)
  
  for(i1 in 1:n_R_se){
    for(i2 in 1:n_beta_se){
      for(i3 in 1:n_lambda_se){
        for(i4 in 1:n_delta_se){
          R = R_se[i1]
          beta = beta_se[i2]
          lambda = lambda_se[i3] 
          delta = delta_se[i4]
          
          index_tem = (i1-1) * n_beta_se * n_lambda_se * n_delta_se
          index_tem = index_tem + (i2-1) * n_lambda_se * n_delta_se
          index_tem = index_tem + (i3-1) * n_delta_se
          index_tem = index_tem + i4
          
          fit = try(vali_middle(R, beta, lambda, delta, xkn_train, xkn_vali, xkn_test,y_train, y_vali, y_test))
          if ("try-error" %in% class(fit)) {
          }
          else {
            res_tem = fit
          }
          y_train_pre = res_tem$y_train_fit
          y_vali_pre = res_tem$y_vali
          y_test_pre = res_tem$y_test
          
          res_mat[index_tem, 1] = sum((y_train_pre - y_train)^2)/n_train
          res_mat[index_tem, 2] = sum((y_vali_pre - y_vali)^2)/n_vali
          res_mat[index_tem, 3] = sum((y_test_pre - y_test)^2)/n_test
          print(index_tem)
          
        }
      }
    }
  }
  

  index = which.min(res_mat[,2])
  MSE_fit = res_mat[index,1]
  MSE_vali = res_mat[index,2]
  MSE_pre = res_mat[index,3]
  res= list(MSE_vali=MSE_vali,MSE_pre=MSE_pre,MSE_fit =MSE_fit, res_mat=res_mat)
  return(res)
}



vali_middle = function(R, beta, lambda, delta, xkn_train, xkn_vali, xkn_test,y_train, y_vali, y_test){
  N = 2000
  K = length(xkn_train[[1]])
  n_train = length(xkn_train)
  n_test = length(xkn_test)
  n_vali = length(xkn_vali)
  
  
  f =  array(rnorm(R*K*n_train), c(R, K, n_train))
  f_test = array(0, c(R,K,n_test))
  f_vali = array(0, c(R,K,n_vali))
  f_train_fit = array(0, c(R,K,n_train))
  
  
  Kmat_all = K_com_all(xkn_train,delta)/lambda
  Kmat_all_inv = Kmat_all
  for(k in 1:K){
    Kmat_all_inv[,,k] = solve(Kmat_all[,,k])
  }
  
  Kmat_test = K_com_pre(xkn_test, xkn_train, delta)/lambda
  Kmat_vali = K_com_pre(xkn_vali, xkn_train, delta)/lambda
  Kmat_train =  K_com_pre(xkn_train, xkn_train, delta)/lambda
  
  y_test = 0
  y_vali = 0
  y_train_fit = 0
  for(iter in 1:N){
    #print(iter)
    for(k in 1:K){
      for(r in 1:R){
        a = c()
        b = c()
        kvec = 1:K
        rvec = 1:R
        for(i in 1:n_train){
          a[i] = prod(f[r,kvec[-k], i])
          tem = 0
          for(r_tem in rvec[-r]){
            tem = tem +  prod(f[r_tem,,i]) 
          }
          b[i] = tem
        }
        Krk = Kmat_all[,,k]
        Sigmark = chol2inv(chol(Kmat_all_inv[,,k] + diag(a^2)/beta))
        murk = Sigmark %*% (hadamard.prod(a,y_train-b))/beta
        f[r,k,] =  mvrnorm(1, murk, Sigmark)
      }
    }
    
    #prediction by posterier mean
    if(iter >=1000){
      
      
      tem_fit = matrix(1, n_train, R)
      for(k in 1:K){
        f_train_fit[,k,] = t(Cprod(Cprod(Kmat_train[,,k], Kmat_all_inv[,,k]),t(matrix(c(f[,k,]),R,n_train))))
        tem_fit = hadamard.prod(tem_fit,t(f_train_fit[,k,])) 
      }
      y_train_fit = y_train_fit + rowSums(tem_fit) 
      
      
      
      tem_vali = matrix(1, n_vali, R)
      for(k in 1:K){
        f_vali[,k,] = t(Cprod(Cprod(Kmat_vali[,,k], Kmat_all_inv[,,k]),t(matrix(c(f[,k,]),R,n_train))))
        tem_vali = hadamard.prod(tem_vali,t(f_vali[,k,])) 
      }
      y_vali = y_vali + rowSums(tem_vali) 
      
      
      
      tem_test = matrix(1, n_test, R)
      for(k in 1:K){
        f_test[,k,] = t(Cprod(Cprod(Kmat_test[,,k], Kmat_all_inv[,,k]),t(matrix(c(f[,k,]),R,n_train))))
        tem_test = hadamard.prod(tem_test,t(f_test[,k,])) 
      }
      y_test = y_test + rowSums(tem_test) 
    }
  }
  y_vali = y_vali/(N-1000)
  y_test = y_test/(N-1000)
  y_train_fit = y_train_fit/((N-1000))
  res= list(y_train_fit=y_train_fit,y_vali = y_vali, y_test=y_test)
  return(res)
}


K_com_all <-function(xkn,delta=15){
  n = length(xkn)
  K = length(xkn[[1]])
  Kmat_all = array(0, c(n, n, K))
  for(k in 1:K){
    for(i in 1:n){
      for(j in 1:i){
        Kmat_all[i,j,k] = - sum((xkn[[i]][[k]] - xkn[[j]][[k]])^2)/(2*delta^2)
        Kmat_all[j,i,k] = Kmat_all[i,j,k] 
      }
    }
  }
  return(exp(Kmat_all))
}

K_com_pre <- function(xkpre, xkn,delta=15){
  n = length(xkn)
  K = length(xkn[[1]])
  n_pre = length(xkpre)
  Kmat_pre = array(0, c(n_pre, n, K))
  for(k in 1:K){
    for(i in 1:n_pre){
      for(j in 1:n){
        Kmat_pre[i,j,k] = - sum((xkpre[[i]][[k]] - xkn[[j]][[k]])^2)/(2*delta^2)
      }
    }
  }
  return(exp(Kmat_pre))
}

#- sum((xi - xj)^2)/(2*15^2)



# 
# n_train = 100
# n_vali = 40
# n_test = 500
# p = c(5,5)
# xkn_train = list()
# xkn_vali = list()
# xkn_test = list()
# 
# X_big_train = array(runif(p[1]*n_train*2), c(p[1],2,n_train))
# #X_big_vali = 10*array(runif(p[1]*n_test*2), c(p[1],2,n_test))
# X_big_test = array(runif(p[1]*n_test*2), c(p[1],2,n_test))
# X_big_vali = X_big_train
# 
# X_train = array(0, c(p,n_train))
# X_test = array(0, c(p,n_test))
# 
# y_train = c()
# y_vali = c()
# y_test = c()
# for(i in 1:n_train){
#   tem = list()
#   tem[[1]] = as.matrix(X_big_train[,1,i])
#   tem[[2]] = as.matrix(X_big_train[,2,i])
#   xkn_train[[i]] = tem
#   y_train[i] = sum(tem[[1]][1:5] * tem[[2]][1:5]) + rnorm(1)
# }
# 
# for(i in 1:n_vali){
#   tem = list()
#   tem[[1]] = as.matrix(X_big_vali[,1,i])
#   tem[[2]] = as.matrix(X_big_vali[,2,i])
#   xkn_vali[[i]] = tem
#   y_vali[i] = sum(tem[[1]][1:5] * tem[[2]][1:5]) + rnorm(1)
# }
# 
# for(i in 1:n_test){
#   tem = list()
#   tem[[1]] = as.matrix(X_big_test[,1,i])
#   tem[[2]] = as.matrix(X_big_test[,2,i])
#   xkn_test[[i]] = tem
#   y_test[i] = sum(tem[[1]][1:5] * tem[[2]][1:5]) + rnorm(1)
# }
#


#ypre= kanagawa_GP(R=2, xkn, xkn, y_train)

# sum((ypre- y_vali)^2)/var(y_vali)
# sum((ypre- y_train)^2)/var(y_train)
# ypre
# sum((ypre - y_test)^2)/var(y_test)

# 
# res =vali_middle(R=2, beta=0.1, lambda=0.1, delta=1, xkn_train, xkn_vali, xkn_test,y_train, y_vali, y_test)
# sum((res$y_train_fit - y_train)^2)/var(y_train)
# sum((res$y_vali - y_vali)^2)/var(y_vali)
# sum((res$y_test - y_test)^2)/var(y_test)

#R_se, beta_se,lambda_se,delta_se


#res= kanagawa_GP_vali(2, c(0.1),0.1,c(0.001,0.1,1,10), xkn_train, xkn_vali, xkn_test,y_train, y_vali, y_test)

#res$res_mat


# D_mat <- function(X_big,X_big_pre){
#   n = dim(X_big)[3]
#   n_pre = dim(X_big_pre)[3]
#   Dmat = matrix(0,n_pre, n)
#   for(i in 1:n_pre){
#     for(j in 1:n){
#       Dmat[i,j] = sum((X_big_pre[,,i] - X_big[,,j])^2)
#     }
#   }
#   # mintem = c()
#   # for(i in 1:n_pre){
#   # mintem[i] = min(Dmqt[i,])
#   # }
#   return(Dmat)
# }
# 
# temtem = D_mat(X_big,X_big_pre)
# temvec = c()
# for(i in 1:100){
#   temvec[i] = min(temtem[i,])
# }
# 
# abs(y_test- ypre)
# temvec
#eigen(K_com_all(xkn,0.01)[,,1])$value
#system.time(chol2inv(matrix(rnorm(25000000),5000,5000)))


# 
# n = 1000
# n_use=500
# n_train=n_use * 0.8
# n_vali=n_use *0.2
# n_test = n - n_use
# id_train=1:n_train
# id_vali=(n_train+1):(n_train+n_vali)
# id_test=501:1000 # test 
# 
# X_data = array(0, c(64,64,n))
# X_kanagawa = array(0, c(64,2,n))
# for(i in 1:n){
#   tem1 = runif(64)
#   tem2 = runif(64)
#   X_data[,,i] = as.matrix(tem1) %*% tem2 
#   X_kanagawa[,1,i] = tem1
#   X_kanagawa[,2,i] = tem2
# }
# 
# xkn_train = list()
# xkn_vali = list()
# xkn_test = list()
# X_kanagawa_train = X_kanagawa[,,id_train]
# X_kanagawa_vali = X_kanagawa[,,id_vali]
# X_kanagawa_test = X_kanagawa[,,id_test]
# for(i in 1:n_train){
#   tem = list()
#   tem[[1]] = as.matrix(X_kanagawa_train[,1,i])
#   tem[[2]] = as.matrix(X_kanagawa_train[,1,i])
#   xkn_train[[i]] = tem
# }
# 
# for(i in 1:n_vali){
#   tem = list()
#   tem[[1]] = as.matrix(X_kanagawa_vali[,1,i])
#   tem[[2]] = as.matrix(X_kanagawa_vali[,1,i])
#   xkn_vali[[i]] = tem
# }
# 
# for(i in 1:n_test){
#   tem = list()
#   tem[[1]] = as.matrix(X_kanagawa_test[,1,i])
#   tem[[2]] = as.matrix(X_kanagawa_test[,1,i])
#   xkn_test[[i]] = tem
# }
# 
# 
# 
# lambda = 10 #在某个例子中，这个越小，拟合的越好
# delta=1
# K=2
# #sampling
# Kmat_all = K_com_all(xkn_train,delta)/lambda
# Kmat_all_inv = Kmat_all
# for(k in 1:K){
#   Kmat_all_inv[,,k] = solve(Kmat_all[,,k])
# }
# 
