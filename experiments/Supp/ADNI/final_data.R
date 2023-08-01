library(png)
library(BNTR)
library(glmnet)
library(rTensor)
library(MASS)
library("RNifti")
source('~/Desktop/ADHD/Athena/Rdata/AddFun.R', echo=TRUE)

#load("~/Desktop/ADNI/dataADNIperson1_80803_51_53.Rdata")
load("~/Desktop/ADNI/data_all80803.Rdata")
sex_age = data_all$sex_age
rmindex= which(sex_age[,1]==0)
X_use = data_all$X_all[-rmindex,]
y_use = data_all$y_all[-rmindex]
sex_age = sex_age[-rmindex,]

n_all = length(y_use)
X_tem = array(0, c(80,80,3,n_all))
for(i in 1:n_all){
  X_tem[,,,i] = array(X_use[i,], c(80,80,3))
}


X_use = X_tem
X_tem = array(0, c(40,40,3,n_all))
for(i in 1:n_all){
  X_tem[,,,i]=resize_array(X_use[,,,i],c(40,40,3))
}

X_use = X_tem
X_ori = X_tem

varX = var_mesure(X_use)

X_use = zero_one_dis(X_use)

load("~/Desktop/ADNI/rmindex.Rdata")
X_all = X_use[,,,-rmindex]
X_ori = X_ori[,,,-rmindex]
y_all = y_use[-rmindex]
sex_age = sex_age[-rmindex,]
n_all = length(y_all)