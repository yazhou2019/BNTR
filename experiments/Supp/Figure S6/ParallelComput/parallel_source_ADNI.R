#parallel sourece
#parallel sourece
library(stringr)
library(rTensor)
library(glmnet)
library(Matrix)
library(MASS)

source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg1220.R')

R=c(1,2,3,4,5,6,7,8)
alpha=c(0,0.5,1)
#alpha=c(0)
lambda=c(0.001,0.0025,0.005,0.0075,0.01,0.025,0.05,0.075,0.1,0.25,0.5,0.75,1,2.5,5,7.5,10,25,50,75,100,250,500,750,1000)




load("/data/yzhou/BNTR/NewReal/data40403.Rdata")


X_all = data$X_all
y_all = data$y_all

rm(data)
gc()

