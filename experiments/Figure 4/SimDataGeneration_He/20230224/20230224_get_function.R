f1_old <-function(x){
  res = x + 0.6 * sin(2*pi*(x-0.5)^2)
  return(res)
}

f2_old <- function(x){
  res = x + 0.3 * cos(2*pi*x)
  return(res)
}


f1 <- function(x){
  res_he = 2*exp(x^2) * x - exp(x)
  res = 0.5*x*exp(x^2) * x - 0.25*exp(x)*x
  return(res)
}

f2 <- function(x){
  res_he = (2*x -1)/(x^2 - x - 2)
  res = (2*x -1)*x/(x^2 - x -2)
  return(res)
}

f3 <- function(x){
  #res = 3*x^2 - 2*x
  res_he = 3*x^2 - 2*x
  res = x^2 - x
  return(res)
}


f4 <- function(x){
  #res = 2^0.5 * sin(2*pi*x)
  res_he = 2^0.5 * sin(2*pi*x)
  res = 0.25 * sin(2*pi*x)
  return(res)
}

f5 <- function(x){
  res_he = sinh(x-0.5)
  res = x*sinh(x-0.5)
  return(res)
}

falls = list()
falls[[1]] = f1
falls[[2]] = f2
falls[[3]] = f3
falls[[4]] = f4
falls[[5]] = f5

plotfunction <-function(f){
  y = c()
  x = c()
  for (i in 1:100){
    x[i] = i/100
    y[i] = f(x[i])
    plot(x, y)
  }
  return(y)
}

int_function <- function(f, precise = 10000){
  y = f(c(0, 1:precise)/precise)
  return(sum(y)/precise)
}

l2norm_function <- function(f){
  precise = 10000
  y = f(c(1:precise)/precise)^2
  return(sqrt(sum(y)/precise))
}

if(1==0){
plotfunction(f1)
plotfunction(f2)
plotfunction(f3)
plotfunction(f4)
plotfunction(f5)
}

print("integral")
print(int_function(f1))
print(int_function(f2))
print(int_function(f3))
print(int_function(f4))
print(int_function(f5))
print("l2norm")
print(l2norm_function(f1))
print(l2norm_function(f2))
print(l2norm_function(f3))
print(l2norm_function(f4))
print(l2norm_function(f5))

#之前正文的函数
print('old integral')
print(int_function(f1_old))
print(int_function(f2_old))
print("old l2norm")
print(l2norm_function(f1_old))
print(l2norm_function(f2_old))

# 
# library(stats)
# integrate(f5,0,1)
# y = plotfunction(f2)
