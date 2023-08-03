load("./data/y_all.Rdata")
y_true = y_all
source('./ComponentsNonLin/functions_needed.R')
load("./data/monkey_1_new_K11_tuning.Rdata")
set.seed(2023)
load("./data/X_test_1.Rdata")
load("./data/X_train_1.Rdata")
X_all=array(0,c(64,10,10,10000))
X_all[,,,1:4000]=X_train
X_all[,,,4001:10000]=X_test
rm(X_train)
rm(X_test)
gc()

num_knots = 10
knots_used = stats::quantile(c(X_all), probs = c(seq(0, 1, 1/(num_knots - 1))))
tildePhiX = tildePhiX_trans(X_all, knots_used)


mse_pre = c()
for (iter in 1:10){
  mse_pre[iter] = result[[iter]][[2]]
}
index = which(order(mse_pre)==6)
BB = full_R(result[[index]][[4]])
y_all = crossprod(matrix(tildePhiX, c(prod(dim(BB)), 10000)), as.vector(BB)) #+ rnorm(10000,0,1)
f_error = y_true -  (y_all + 45)

set.seed(432)
y_new = y_all + f_error[sample(1:10000,10000, replace = T)]
#print(var(y_new))

y_all = y_new 
setwd("./data")
save(y_all,file = "y_all_sim_new2.Rdata")
#y_all_sim_new.Rdata对应 y_new = y_all + f_error[sample(1:10000,10000)],set.seed(2022)
setwd("../")


# generate_sim_real_with_norm_noise <- function(){
#   setwd("~/Desktop/Research/JRSSB/upload2_newserver/BNTR/RealDataSim")
#   load("./data/y_all.Rdata")
#   y_true = y_all
#   source('./ComponentsNonLin/functions_needed.R')
#   load("./data/monkey_1_new_K11_tuning.Rdata")
#   
#   set.seed(2023)
#   load("./data/X_test_1.Rdata")
#   load("./data/X_train_1.Rdata")
#   X_all=array(0,c(64,10,10,10000))
#   X_all[,,,1:4000]=X_train
#   X_all[,,,4001:10000]=X_test
#   rm(X_train)
#   rm(X_test)
#   gc()
#   
#   num_knots = 10
#   knots_used = stats::quantile(c(X_all), probs = c(seq(0, 1, 1/(num_knots - 1))))
#   tildePhiX = tildePhiX_trans(X_all, knots_used)
#   
#   
#   mse_pre = c()
#   for (iter in 1:10){
#     mse_pre[iter] = result[[iter]][[2]]
#   }
#   index = which(order(mse_pre)==6)
#   BB = full_R(result[[index]][[4]])
#   y_all = crossprod(matrix(tildePhiX, c(prod(dim(BB)), 10000)), as.vector(BB)) #+ rnorm(10000,0,1)
#   
#   set.seed(432)
#   y_new = y_all + rnorm(10000, mean=0, sd=1.2*sqrt(var(y_all)))
#   #f_error = y_true -  (y_all + 45)
#   
#   #set.seed(432)
#   #y_new = y_all + f_error[sample(1:10000,10000, replace = T)]
#   print(var(y_new))
#   
#   y_all = y_new 
#   setwd("./data")
#   #save(y_all,file = "y_all_sim_new3.Rdata")
#   save(y_all,file = "y_all_sim_new_SNR120.Rdata")
#   #y_all_sim_new.Rdata对应 y_new = y_all + f_error[sample(1:10000,10000)],set.seed(2022)
#   setwd("../")
#   
# }












