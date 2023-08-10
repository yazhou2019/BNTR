library(stringr)
library(ggplot2)
library(cowplot)
library(showtext)





n_use_se = c(500,1000)

ppp=list()
for(i in 1:2){
n_use = n_use_se[i]
name = str_c("~/Desktop/JRSSBreview/smaller objective value/res", n_use, "_smaller_oj.Rdata")
load(name)

dataplot = data.frame(x=c(rep("LHSloss",50), rep("RHSobj",50)), 
                      y=c(log(res$res_mat[,2]), log(res$res_mat_app[,2]))
                      )
  
pp=ggplot(data = dataplot) + geom_boxplot(aes(x = x, y = y)) #+   geom_point(aes(x = x, y = y,color=x))
pp = pp + labs(x=str_c('n=',n_use), y="Log Values",title = str_c("Small Loss Value")) + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(vjust=0.1))
ppp[[i]] = pp
}



for(i in 1:2){
  n_use = n_use_se[i]
  name = str_c("~/Desktop/JRSSBreview/smaller objective value/res", n_use, "_smaller_oj.Rdata")
  load(name)
  
  dataplot = data.frame(x=c(rep("RHSobj-LHSloss",50)), 
                        y=c(log(res$res_mat_app[,2]-res$res_mat[,2]))
  )
  
  pp=ggplot(data = dataplot) + geom_boxplot(aes(x = x, y = y)) #+   geom_point(aes(x = x, y = y,color=x))
  pp = pp + labs(x=str_c('n=',n_use), y="Log Values",title = str_c("Value difference")) + theme(plot.title = element_text(hjust = 0.5),axis.text.x = element_text(vjust=0.1))
  pp = pp + scale_y_continuous(limits=c(0,14))#breaks = c(0, 4, 8, 12, 16))
  ppp[[i+2]] = pp
}


#gg= ggdraw() +  draw_plot(ppp[[1]], 0, 0, 0.5, 1)+draw_plot(ppp[[2]], 0.5, 0, 0.5, 1) #+ draw_plot(ppp[[2]], 0.5, 0, 0.5, 0.5)

gg= ggdraw() +  draw_plot(ppp[[1]], 0, 0.5, 0.5, 0.5)+draw_plot(ppp[[2]], 0.5, 0.5, 0.5, 0.5) + draw_plot(ppp[[3]], 0, 0, 0.5, 0.5) + draw_plot(ppp[[4]], 0.5, 0, 0.5, 0.5)


showtext_begin()
print(gg)
showtext_end()

# png("smaller_loss_500.png",units="in", width=4.5*2, height=4.5*4.03/5.83,res=600)
# gg
# dev.off()

png("smaller_loss_500_new.png",units="in", width=4.5*2, height=4.5*4.03/5.83*2,res=600)
gg
dev.off()
