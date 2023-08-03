
# TLR-rescaled 
mse_mat = matrix(0, 10, 2)
for(iter in 1:10){
  mse_mat[iter,1:2] = result[[iter]][[2]]
}
print(mean(mse_mat[,2]))
print(var(mse_mat[,2]))

#BNTR
mse_mat = matrix(0, 10, 2)
for(iter in 1:10){
  mse_mat[iter,1:2] = c(result[[iter]][[1]],result[[iter]][[2]])
}
print(mean(mse_mat[,2]))
print(var(mse_mat[,2]))

#ENet
mse_mat = matrix(0, 10, 2)
for(iter in 1:10){
  mse_mat[iter,1:2] = c(result[[iter]][[iter]][[11]],result[[iter]][[iter]][[11]])
}
print(mean(mse_mat[,2]))
print(var(mse_mat[,2]))


