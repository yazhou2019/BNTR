
library(MASS)
library(matrixcalc)
library(BNTR)
# Sigma <- matrix(c(10,3,3,2),2,2)
# Sigma

# mydata <- mvrnorm(n=500,rep(0, 2),Sigma)


kanagawa_GP <- function(R=1, xkn, xkpre, y){
  N = 200
  n = length(xkn)
  K = length(xkn[[1]])
  n_pre = length(xkpre)
  beta = 1
f =  array(rnorm(R*K*n), c(R, K, n))
fpre = array(0, c(R,K,n_pre))
lambda = 10
delta=0.1
#sampling
# f[r,k,] =  mvrnorm(1,rep(0, 2),Sigma)


Kmat_all = K_com_all(xkn,delta)/lambda
Kmat_all_inv = Kmat_all
for(k in 1:K){
  #Kmat_all_inv[,,k] = chol2inv(Kmat_all[,,k])
  Kmat_all_inv[,,k] = solve(Kmat_all[,,k])
}
Kmat_pre = K_com_pre(xkpre, xkn, delta)/lambda

ypre = 0
for(iter in 1:N){
for(k in 1:K){
  for(r in 1:R){
#caculation a,b
#1, 先计算a, b
a = c()
b = c()
kvec = 1:K
rvec = 1:R
for(i in 1:n){
  a[i] = prod(f[r,kvec[-k], i])
  tem = 0
  for(r_tem in rvec[-r]){
  tem = tem +  prod(f[r_tem,,i]) 
  }
  b[i] = tem
}

#因为对每个mode都是同一个kernal

Krk = Kmat_all[,,k]
#print(a)
#2，再计算协方差矩阵和均值
print(a)
#Sigmark = Krk - Krk %*% chol2inv(Krk + diag(a^2)/beta) %*% Krk
#Sigmark = Krk - Krk %*% ginv(Krk + diag(a^2)/beta) %*% Krk
Sigmark = chol2inv(Kmat_all_inv[,,k] + beta *diag(1/a^2))

murk = Sigmark %*% (hadamard.prod(a,y-b))/beta
#print(Sigmark)
#print(hadamard.prod(a,y-b))
#3 抽样，这里加0.00001*diag(n)是为了预防协方差矩阵为负
#f[r,k,] =  mvrnorm(1, murk, Sigmark+0.00001*diag(n))
#eigen(Sigmark)$value
#print(eigen(Sigmark)$value)
f[r,k,] =  mvrnorm(1, murk, Sigmark +0.00001*diag(n))
}
}

tem = matrix(1, n_pre, R)
for(k in 1:K){
  fpre[,k,] = t(Kmat_pre[,,k] %*% Kmat_all_inv[,,k] %*% t(matrix(c(f[,k,]),R,n)))
  tem = tem * t(fpre[,k,]) #n_pre \times R
}
#print(tem)
 ypre = ypre + rowSums(tem) 
 #print(ypre)
  
}
ypre = ypre/N
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
p = c(2,2)
xkn = list()
xkn_pre = list()
X_big = array(rnorm(p[1]*n_train*2), c(p[1],2,n_train))
X_big_pre = array(rnorm(p[1]*n_test*2), c(p[1],2,n_test))

X_train = array(0, c(p,n_train))
X_test = array(0, c(p,n_test))
y = c()
for(i in 1:n_train){
  tem = list()
  tem[[1]] = as.matrix(X_big[,1,i])
  tem[[2]] = as.matrix(X_big[,2,i])
  xkn[[i]] = tem
  X_train[,,i] = full_R(tem)
  y[i] = sum(tem[[1]]) * sum(tem[[2]])
}
y_test =c()
for(i in 1:n_test){
  tem = list()
  tem[[1]] = as.matrix(X_big_pre[,1,i])
  tem[[2]] = as.matrix(X_big_pre[,2,i])
  xkn_pre[[i]] = tem
  X_test[,,i] = full_R(tem)
  y_test[i] = sum(tem[[1]]) * sum(tem[[2]])
}

ypre= kanagawa_GP(R=3, xkn, xkn_pre, y)
ypre
sum((ypre - y_test)^2)/var(y_test)

#eigen(K_com_all(xkn,0.01)[,,1])$value
#system.time(chol2inv(matrix(rnorm(25000000),5000,5000)))