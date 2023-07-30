


#Plot for comparison
library(BNTR)
##All the plot
p=c(64,64)
par(mar = c(0,0,0,0),mfcol=c(1,1),mai=c(0.25,0.3,1,0.25))


color_used = c(5/255,120/255,5/255)
BB = matrix(0, p[1], p[2])

BBplot_bar = array(1, c(p[1],p[2],3))
for(i in 1:p[2]){
  BB[,i]=(i-1)/63
}
BB = BB*2.5
BB[BB>1]=1

BBplot_bar[,,1] = color_used[1]  * ((-1/color_used[1] + 1) *BB + 1/color_used[1])
BBplot_bar[,,2] = color_used[2]  * ((-1/color_used[2] + 1) *BB + 1/color_used[2]) 
BBplot_bar[,,3] = color_used[3] * ((-1/color_used[3] + 1) *BB + 1/color_used[3])





plot(c(3,p[1]-2),c(3,p[2]-2),xaxt="n",yaxt='n',type="n")

library(stringr)
space_use = " "
for(i in 1:25){
  space_use = str_c(space_use, " ")
}
index_use = "0"
for(i in 2:5){
  index_use = str_c(index_use, space_use, (i-1) * 0.25)
}

mtext(index_use,3,line=0.2,cex=3)
#par(cex=5)

#mtext("0                            0.25                            0.5                            0.75                            1 ",3,line=0.2)
rasterImage(BBplot_bar, 1-0.5, 1-0.5, p[1]+1, p[2]+1)
