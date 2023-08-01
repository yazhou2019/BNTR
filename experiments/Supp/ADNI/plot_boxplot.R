library("RNifti")
library(BNTR)
source('~/Desktop/ADHD/Athena/Rdata/AddFun.R', echo=TRUE)
library(rTensor)
library(refund.wave)
library(refund)
library(fields)
library(stringr)


#see the result
load("/Users/zhouya/Desktop/resADNI/resall.Rdata")
mis_error = matrix(0,10, 2)
tuning_par = matrix(0, 10, 6)
for(i in 1:10){
  mis_error[i,1] =resall$TLR2[[i]][[1]][2]
  mis_error[i,2] =resall$Nonlin[[i]][[1]][2]
  tem = resall$TLR2[[i]][[2]]
  index = which.max(tem[,2])
  tuning_par[i,1:2] = tem[index,4:5] 
  tuning_par[i,3] = tem[index,6] * 495
  tem = resall$Nonlin[[i]][[2]]
  BBuse = resall$Nonlin[[i]][[3]]
  index = which.max(tem[,2])
  tuning_par[i,4:6] = tem[index,4:6]
}

library(stringr)
library(ggplot2)
library(cowplot)
library(showtext)

dataplot =data.frame(x=c(rep("TLR_rescaled", 10), rep("BroadcasTR", 10)), y=c(mis_error[,1],mis_error[,2]))
pp=ggplot(data = dataplot) + geom_boxplot(aes(x = x, y = y)) #+   geom_point(aes(x = x, y = y,color=x))
pp = pp + labs(x='', y="Accuracy") + theme(plot.title = element_text(hjust = 0.5))
pp
  

png("TLR_BNTR_com_ADNI.png",units="in", width=4.5, height=4.5*4.03/5.83,res=600)
pp
dev.off()
  