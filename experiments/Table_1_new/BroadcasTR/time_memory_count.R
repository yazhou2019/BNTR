library(stringi)
library(stringr)


setwd("~/Desktop/Research/JRSSB/换editor上传/Table_1")



methods = c("TLR1", "TLR-rescaled","ENetR", "BroadcasTR")
methods_prime = c("TLR1","TLR2", "ENet", "BNTR")
#options(digits=3) 


for(n in c(500,750,1000)){

i=4
method = methods[i]
method_prime = methods_prime[i]
pathname = str_c("./", method, "/SimResults/",method_prime,"_", n, "_20230918.Rdata" )
load(pathname)

time_seq = c()
mem_seq = c()
max_seq = c()
for(iter in 1:50){
  time_seq[iter] =  result[[iter]][[31]]
  mem_seq[iter] = result[[iter]][[32]]
  max_seq[iter] = result[[iter]][[33]]
}
out = str_c("n=", n, " ", round(mean(time_seq),3), " (", round(sqrt(var(time_seq)),3), ")")
print(out)
}

