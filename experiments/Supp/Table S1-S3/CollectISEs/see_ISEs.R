

n_all=c(500,750,1000)
method_all=c("TLR1", "TLR2", "ENet", "BNTR")

n=n_all[3]


hehe_matrix=matrix(0,4,2*4)
for(i in 1:4){
method=method_all[i]

path=str_c("final_error_", method, "_", n, "_new.Rdata")
load(path)


for(signal_i in 1:4){
  hehe_matrix[signal_i,2*(i-1)+1]=mean(final_error[signal_i,])
  hehe_matrix[signal_i,2*(i-1)+2]=sqrt(var(final_error[signal_i,]))
}

}
round(hehe_matrix, 3)

(hehe_matrix[,1]-hehe_matrix[,7])/hehe_matrix[,1]
(hehe_matrix[,3]-hehe_matrix[,7])/hehe_matrix[,3]
(hehe_matrix[,5]-hehe_matrix[,7])/hehe_matrix[,5]
