#true_signal horse, cross, cube
library(stringr)
load('./SimDataGeneration_He/shape/X_cross.Rdata')
load('./SimDataGeneration_He/shape/X_horse.Rdata')
load('./SimDataGeneration_He/shape/X_cube.Rdata')
load('./SimDataGeneration_He/shape/X_mouth.Rdata')
load('./SimDataGeneration_He/shape/X_H.Rdata')
load('./SimDataGeneration_He/shape/butterfly.Rdata')
source('./SimDataGeneration_He/20230224_get_function.R')


get_BB <-function(){

p=c(64,64)
p_small=dim(X_horse)

p1start=(p[1]-p_small[1])/2
p1end=p1start+p_small[1]-1
p2start=(p[2]-p_small[2])/2
p2end=p2start+p_small[2]-1
#X_horse[13:14,6:19]=1
#X_horse[16:18,1:5]=1
#X_horse[12:14,6:8]=1
X_horse[,3:5]=0
X_cross_aug1=matrix(0,p[1],p[2])
X_cross_aug1[p1start:p1end,p2start:p2end]=X_cross

#X_cross_aug2=matrix(0,p[1],p[2])
#X_cross_aug2[1:p_small[1],1:p_small[2]]=X_cross

X_cross_aug2=matrix(0,p[1],p[2])
X_cross_aug2[(p[1]-p_small[1]+1):p[1],(p[2]-p_small[2]+1):p[2]]=X_cross

X_H_aug1=matrix(0,p[1],p[2])
X_H_aug1[1:p_small[1],1:p_small[2]]=X_H

X_mouth_aug1=matrix(0,p[1],p[2])
X_mouth_aug1[1:p_small[1],1:p_small[2]]=X_mouth

X_horse_aug1=matrix(0,p[1],p[2])
X_horse_aug1[p1start:p1end,p2start:p2end]=X_horse

X_cube_aug1=matrix(0,p[1],p[2])
X_cube_aug1[(p[1]-p_small[1]+1):p[1],(p[2]-p_small[2]+1):p[2]]=X_cube




BB_signal_all=list()
BB_signal_all[[1]]<-X_cross_aug1
BB_signal_all[[2]]<-X_horse_aug1
BB_signal_all[[3]]<-X_mouth_aug1 
BB_signal_all[[4]]<-X_cross_aug2
BB_signal_all[[5]]<-X_H_aug1
BB_signal_all[[6]]<-X_cube_aug1
BB_signal_all[[7]]<-butterfly
BB_signal_all[[8]]=array(0, c(64,64))
BB_signal_all[[8]][3:10, 3:10] = matrix(rnorm(8*8, 0, 0.1))
BB_signal_all[[9]]=array(0, c(64,64))
BB_signal_all[[9]][55:62, 55:62] = matrix(rnorm(8*8, 0, 0.1))

#BB_signal_all[[2]] = BB_signal_all[[2]] + BB_signal_all[[7]] + BB_signal_all[[8]]
#100,107
set.seed(100)
a = runif(64, 0.5, 1)
b = runif(64, 0.5, 1)
c = runif(64, -0.5, 0.5)
c[c>0]=1
c[c<=0]=-1
d = runif(64, -0.5, 0.5)
d[d>0]=1
d[d<=0]=-1
a = a*c
b = b*d
mask_mat = matrix(a,64, 1) %*% matrix(b, 1, 64)

for(i in 1:7){
  BB_signal_all[[i]] = BB_signal_all[[i]] * mask_mat
}


BBprod=list()
for(i in 1:9){
  BBprod[[i]]=array(0,c(1,dim(BB_signal_all[[i]])))
  BBprod[[i]][1,,]=BB_signal_all[[i]]
}


BB=list()

BB[[1]]=X_cross_aug1
BB[[2]]=X_horse_aug1
BB[[3]]=X_mouth_aug1+X_cross_aug2
BB[[4]]=X_H_aug1+X_cube_aug1
return(BBprod)
}

get_true_l2norm <- function(BBprod){
  
f1l2norm = (l2norm_function(f1))
f2l2norm = (l2norm_function(f2))
f3l2norm = (l2norm_function(f3))
f4l2norm = (l2norm_function(f4))
f5l2norm = (l2norm_function(f5))

BBprod_l2norm = list()
BBprod_l2norm[[1]] = 0.5 * BBprod[[1]]
BBprod_l2norm[[2]] = BBprod[[2]] * f2l2norm + BBprod[[8]]* f4l2norm  + BBprod[[9]] * f5l2norm
BBprod_l2norm[[3]] = BBprod[[3]] * f1l2norm + BBprod[[4]] * f1l2norm
BBprod_l2norm[[4]] = BBprod[[5]] * f1l2norm + BBprod[[6]] * f3l2norm
BBprod_l2norm[[5]] = BBprod[[7]] * f3l2norm



for(i in 1:5){
  BBprod_l2norm[[i]] = abs(BBprod_l2norm[[i]])  
}

# get_outputs <- function(v, BBprod, X_data){
#   y_all = list()
#   y_all[[1]]=v[[1]]+t(ctprod(BBprod[[1]],X_data,2))
#   y_all[[2]]=v[[2]]+t(ctprod(BBprod[[2]],f1(X_data),2))+t(ctprod(BBprod[[8]],f2(X_data),2))+t(ctprod(BBprod[[9]],f5(X_data),2))
#   y_all[[3]]=v[[3]]+t(ctprod(BBprod[[3]],f2(X_data),2))+t(ctprod(BBprod[[4]],f2(X_data),2))
#   y_all[[4]]=v[[4]]+t(ctprod(BBprod[[5]],f3(X_data),2))+t(ctprod(BBprod[[6]],f5(X_data),2))
#   y_all[[5]]=v[[5]]+t(ctprod(BBprod[[7]],f4(X_data),2))
#   
#   return(y_all)
#}


return(BBprod_l2norm)


}
