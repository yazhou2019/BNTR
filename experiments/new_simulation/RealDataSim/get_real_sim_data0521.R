source('./ComponentsNonLin/functions_needed.R') # import the dependency
load("./data/y_all.Rdata") #  load real observations (outputs, y)
load("./data/fitted_BB.Rdata") # load the fitted BB
load("./data/fitted_b0.Rdata") # load the fitted B0
load("./data/X_test_1.Rdata") # load real observations (inputs, X)
load("./data/X_train_1.Rdata") # load real observations (inputs, X)
load("./data/idtest_matrix.Rdata") # random splits 
load("./data/idtrain_matrix.Rdata") # random splits

# X_all are the real observations 
X_all=array(0,c(64,10,10,10000)) 
X_all[,,,1:4000]=X_train
X_all[,,,4001:10000]=X_test
rm(X_train)
rm(X_test)
gc() # save the memory

# generate the truncated power basis
knots_used = stats::quantile(c(X_all), probs = c(seq(0, 1, 1/(10 - 1))))
tildePhiX = tildePhiX_trans(X_all, knots_used)

# obtain the error and hat(m)
y_sim = crossprod(matrix(tildePhiX, c(prod(dim(BB)), 10000)), as.vector(BB))  + b0 
f_error = y_all - y_sim 


# obtain ten data splits
get_splits <- function(idx){
  # idx = 1,2,3,4,5,6,7,8,9,10; corresponding to 10 replications
  id_train = idtrain_matrix[idx,]
  id_vali = idtest_matrix[idx, 1:1000]
  id_test = idtest_matrix[idx, 1001:6000]
  id_list = list(id_train=id_train, id_test=id_test, id_vali=id_vali)
}
rep_splits = list()
for(idx in 1:10){
  rep_splits = get_splits(idx)
}


set.seed(432)
# generate the simulated response
y_sim = y_sim + f_error[sample(1:10000,10000, replace = T)]


##############
##############
#rep_splits: ten random splits
#y_sim: the simulated response

