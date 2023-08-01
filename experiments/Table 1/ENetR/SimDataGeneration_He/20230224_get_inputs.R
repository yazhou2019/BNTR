


library(mnormt)
get_inputs <- function(n=10, p = c(32,32), distri="sen_norm"){
  
  if(distri=="sen_norm"){
    varcov = diag(p[1] * p[2])
    rho = 0.5
    for (i in 1:p[1]^2){
      for (j in 1:p[2]^2){
        i1 = i%/%p[1]
        j1 = i%%p[1]
        i2 = j%/%p[2]
        j2 = j%%p[2]
        varcov[i,j] = rho^(abs(i1-i2)+abs(j1-j2))
      }
    }
    
    X =  rmnorm(n = n, mean = rep(0.5,p[1]*p[2]), varcov=varcov)
    X[X>1]=1
    X[X<0]=0
    X_tem = array(0, c(p, n))
    for (i in 1:n){
      X_tem[,,i] = matrix(X[i,], p)
    }
    X = X_tem 
  }else if(distri=="uniform"){
    
    X = array(runif(n * p[1] * p[2], 0, 1), c(p, n))
  }
  
  return(X)
}
#X = get_inputs(n=10)