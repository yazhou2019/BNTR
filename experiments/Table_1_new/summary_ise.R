library(stringi)
library(stringr)


setwd("~/Desktop/Research/JRSSB/换editor上传/Table_1")


n = 500
methods = c("TLR1", "TLR-rescaled","ENetR", "BroadcasTR")
methods_prime = c("TLR1","TLR2", "ENet", "BNTR")
options(digits=3) 

for (case in 1:5){
res = c(" ")
for(i in 1:4){
method = methods[i]
method_prime = methods_prime[i]
pathname = str_c("./", method, "/SimResults/ise_",method_prime,"_", n, "_20230918.Rdata" )
load(pathname)
seq = ise_mat[,2 * case - 1]
res = str_c(res, "&", round(mean(seq),3)," (", round(sqrt(var(seq)),3), ") ")
if(i==4){
  res = str_c(res, "\\")
}
}
print(res)
}
