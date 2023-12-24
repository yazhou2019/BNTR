#parallel computing for multiple replications


library(snowfall)
sfInit(parallel=TRUE, cpus=25, type="SOCK")




wrapper <- function(idx) {
  # Output progress in worker logfile
  cat( "Current index: ", idx, "\n" )
  #all the data set
  set.seed(idx*100)
  n=1000
  
  ##########################user define in the parallel_source########################################
  n_use=n_use
  ##########################user define in the parallel_source########################################
  n_train=n_use * 0.8
  n_vali=n_use *0.2
  
  n_all = n_train + n_vali
  
  #indicator = 0
  #sigma = 0
  #distri = "sen_norm" #"sen_norm" #"uniform" #
  #p_test = c(64, 64)
  #X_data = get_inputs(n=n_all, p=p_test, distri = distri)
  
  #y_all=list()
  #for (CASE in 1:5){
  #  BB = get_BB_together(case=CASE)
  #  data_list = get_outputs(X=X_data, case=CASE, indicator =  indicator)
  #  y_all[[CASE]] = data_list$y +  get_rd_error(n=n_all, case=CASE, distri = distri, indicator =  indicator, sigma=sigma)
  #}
  
  #load(str_c("./SimResults/BNTR", n_use, "_new_K8_20230302.Rdata"))
  load(str_c("./SimResults/BNTR", n_use, "_20230918.Rdata"))
  #load(str_c("./SimResults/Simul_n", n_use,,"_rep50_final_fix_new.Rdata"))
  
  #load("./SimResults/BNTR1000_new_K8_20230302.Rdata")
  #load("./SimResults/Simul_n1000_rep50_final_fix_new.Rdata")
  load("./SimResults/Simul_n1000_rep50_final_fix_new_editor.Rdata")
  
  BB_list = data_all[[51]]
  
  #num_knots=7
  
  ise_mat = matrix(0, 2, 5)
  for(signal_i in 1:5){
    #id_train=1:n_train
    #id_vali=(1000-n_vali+1):1000
    #id_test=851:1000 # test is useless here
    
    #knots_used = stats::quantile(c(X_data[,,id_train]), probs = c(seq(0, 1, 1/(num_knots - 1))))
    
    BB = result[[idx]][[signal_i]]
    bb = result[[idx]][[signal_i+5]]
    MSE_pre = result[[idx]][[signal_i+10]]
    knots_used = result[[idx]][[signal_i+15]]
    BBv = result[[idx]][[signal_i+20]]
    vv = result[[idx]][[signal_i+25]] #[[signal_i]]
    ise_mat[1, signal_i] = final_ise(b0=bb, BB=BB, BBv= BB_list, v =vv, case=signal_i, precision=1000, knots = knots_used, order = 4)$ise
    ise_mat[2, signal_i] = MSE_pre
    
  }
  
  
  return(ise_mat)
}

sfSource('./ParallelComput/parallel_source_500.R')

sfClusterSetupRNG()

result <- sfLapply(1:50, wrapper)

ise_mat= matrix(0, 50, 10)
for(iter in 1:50){
  ise_mat[iter, c(1,3,5,7,9)] = result[[iter]][1,]
  ise_mat[iter, c(2,4,6,8,10)] = result[[iter]][2,]
}

sfStop()
filename <- str_c("ise_BNTR_", n_use, "_20230918.Rdata")
setwd("./SimResults")
save(ise_mat,file = filename)
setwd("../")
