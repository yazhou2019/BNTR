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
library(BNTR)
# Sigma <- matrix(c(10,3,3,2),2,2)
# Sigma

# mydata <- mvrnorm(n=500,rep(0, 2),Sigma)


kanagawa_GP <- function(R=1, xkn, xkpre, y){
  N = 5000
  n = length(xkn)
  K = length(xkn[[1]])
  n_pre = length(xkpre)
  beta = 1000 #在某个例子中，这个越小，拟合的越好
  f =  array(rnorm(R*K*n), c(R, K, n))
  fpre = array(0, c(R,K,n_pre))
  
  lambda = 0.001
  delta=1
  #sampling
  # f[r,k,] =  mvrnorm(1,rep(0, 2),Sigma)
  
  
  Kmat_all = K_com_all(xkn,delta)/lambda
  Kmat_all_inv = Kmat_all
  for(k in 1:K){
    #Kmat_all_inv[,,k] = chol2inv(Kmat_all[,,k])
    #Kmat_all_inv[,,k] = solve(Kmat_all[,,k])
    Kmat_all_inv[,,k] = solve(Kmat_all[,,k])
  }
  
  Kmat_pre = K_com_pre(xkpre, xkn, delta)/lambda
  
  ypre = 0
  for(iter in 1:N){
    print(iter)
    for(k in 1:K){
      for(r in 1:R){
        #caculation a,b
        #1, 先计算a, b
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
        
        #因为对每个mode都是同一个kernal
        Krk = Kmat_all[,,k]
        #print(a)
        #2，再计算协方差矩阵和均值
       # Sigmark = chol2inv(chol(Kmat_all_inv[,,k] + beta *diag(1/a^2)))
        Sigmark = chol2inv(chol(Kmat_all_inv[,,k] + diag(a^2)/beta))
        
        #Sigmark = Krk - Cprod(Cprod(Krk,chol2inv(Krk + diag(a^2)/beta)), Krk)
        #Sigmark = Kmat_all_inv[,,k] 
        #Sigmark = diag(n)
        
        murk = Sigmark %*% (hadamard.prod(a,y-b))/beta
        #print(hadamard.prod(a,y-b))
        #print(murk)
        f[r,k,] =  mvrnorm(1, murk, Sigmark)
      }
    }
    
    if(iter >=1000){
    tem = matrix(1, n_pre, R)
    for(k in 1:K){
      #f[,k,]的维度是R \times n
      #fre[,k,] 的维度是 R\times n
      fpre[,k,] = t(Kmat_pre[,,k] %*% Kmat_all_inv[,,k] %*% t(matrix(c(f[,k,]),R,n)))
      tem = hadamard.prod(tem,t(fpre[,k,])) #n_pre \times R
    }
    #print(tem)
    ypre = ypre + rowSums(tem) 
    #print(ypre)
    
    }
  }
  ypre = ypre/(N-1000)
  return(ypre)
}
#input 形式


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


n_train = 400
n_test = 100
p = c(10,10)
xkn = list()
xkn_pre = list()

X_big = 10*array(runif(p[1]*n_train*2), c(p[1],2,n_train))
X_big_pre = 10*array(runif(p[1]*n_test*2), c(p[1],2,n_test))

X_train = array(0, c(p,n_train))
X_test = array(0, c(p,n_test))
y_train = c()
for(i in 1:n_train){
  tem = list()
  tem[[1]] = as.matrix(X_big[,1,i])
  tem[[2]] = as.matrix(X_big[,2,i])
  xkn[[i]] = tem
  X_train[,,i] = full_R(tem)
  y_train[i] = tem[[1]][1] * tem[[2]][1] + tem[[1]][2] * tem[[2]][2]
}
y_test =c()
for(i in 1:n_test){
  tem = list()
  tem[[1]] = as.matrix(X_big_pre[,1,i])
  tem[[2]] = as.matrix(X_big_pre[,2,i])
  xkn_pre[[i]] = tem
  X_test[,,i] = full_R(tem)
  y_test[i] = tem[[1]][1] * tem[[2]][1] + tem[[1]][2] * tem[[2]][2]
}

ypre= kanagawa_GP(R=2, xkn, xkn_pre, y_train)

sum((ypre- y_train)^2)/var(y_train)
ypre
sum((ypre - y_test)^2)/var(y_test)

#eigen(K_com_all(xkn,0.01)[,,1])$value
#system.time(chol2inv(matrix(rnorm(25000000),5000,5000)))