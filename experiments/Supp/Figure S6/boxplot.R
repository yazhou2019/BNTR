
load("./RealResults/TLR2_linear_ADNI40403.Rdata")
TLR2 = matrix(0, 10, 2)
for(i in 1:10){
  TLR2[i,] = result[[i]][[1]]
}

load("./RealResults/ADNI40403.Rdata")
Nonlin = matrix(0, 10, 2)
for(i in 1:10){
  Nonlin[i,] = result[[i]][[1]]
}

resall = list(TLR2 = TLR2, Nonlin =Nonlin)
save(file="resall.Rdata", resall)

mis_error = matrix(0,10, 2)
mis_error[,1] = resall$TLR2[,2]
mis_error[,2] = resall$Nonlin[,2]

library(stringr)
library(ggplot2)
library(cowplot)
library(showtext)

dataplot =data.frame(x=c(rep("TLR_rescaled", 10), rep("BroadcasTR", 10)), y=c(mis_error[,1],mis_error[,2]))
pp=ggplot(data = dataplot) + geom_boxplot(aes(x = x, y = y)) #+   geom_point(aes(x = x, y = y,color=x))
pp = pp + labs(x='', y="Accuracy") + theme(plot.title = element_text(hjust = 0.5))
pp

png("TLR2NonLinBoxplot.png",units="in", width=4.5, height=4.5*4.03/5.83,res=600)
pp
dev.off()









