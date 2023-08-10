#true_signal horse, cross, cube
library(stringr)
load('./SimDataGeneration/X_cross.Rdata')
load('./SimDataGeneration/X_horse.Rdata')
load('./SimDataGeneration/X_cube.Rdata')
load('./SimDataGeneration/X_mouth.Rdata')
load('./SimDataGeneration/X_H.Rdata')
load('./SimDataGeneration/butterfly.Rdata')


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


BB=list()

BB[[1]]=X_cross_aug1
BB[[2]]=X_horse_aug1
BB[[3]]=X_mouth_aug1+X_cross_aug2
BB[[4]]=X_H_aug1+X_cube_aug1
BB[[5]] = butterfly

# see the true signal
# par(mar = c(0,0,0,0),mfrow=c(3,4),mai=c(0.25,0.25,0.25,0.25))

# for(signal_i in 1:length(BB)){
#  title=str_c('True signal ',signal_i)
#
# plot(c(1,p[1]),c(1,p[2]),xaxt="n",yaxt='n',type="n")
#  axis(side=1,at=c(1,round(p[1]/2),p[1]))
#  axis(side=2,at=c(1,round(p[2]/2),p[2]))
#  rasterImage(BB[[signal_i]], 1, 1, p[1], p[2])
#  title(main =title)
#}

