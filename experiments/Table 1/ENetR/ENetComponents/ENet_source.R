library(stringr)
library(stringi)
source('./ENetComponents/Validation_Enet.R')
source('./SimDataGeneration_He/20230224_get_ise.R')
alpha=c(0,0.5,1)
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)

n=1000
n_use=750
n_train=n_use*0.8
n_vali=n_use*0.2
