library(snowfall)
sfInit(parallel=TRUE, cpus=25, type="SOCK")



wrapper <- function(idx) {
  # Output progress in worker logfile
  cat( "Current index: ", idx, "\n" )
  #BB_idx=list()
  #all the data set
  n_seq_list = list()
  n_seq = c(250, 500, 750, 1000)

  for(samplesize_i in 1:4){
  n=n_seq[samplesize_i]
  n_data=n
  p=c(64,64)
  if(n_data==1000){
    num_knots=6
  }
  if(n_data==750||n_data==500){
    num_knots=5
  }
  if(n_data==250){
    num_knots=4  
  }  
  if(n==250){
    path1 = "/data/yzhou/BNTR/FullSimPaper/SimResults/BNTR250_tnorm_K4.Rdata"
  }
  if(n==500){
    path1 = "/data/yzhou/BNTR/FullSimPaper/SimResults/BNTR500_tnorm_K5.Rdata"
  }
  if(n==750){
    path1 = "/data/yzhou/BNTR/FullSimPaper/SimResults/BNTR750_tnorm_K5.Rdata"
  }
  if(n==1000){
    path1 = "/data/yzhou/BNTR/FullSimPaper/SimResults/BNTR1000_tnorm_K7.Rdata"
  }
  
  load(path1)
  set.seed(idx*100)
  X = array(rtnorm(64*64*1000, mean=0.5, sd=0.25, lower=0, upper=1),c(64,64,1000))
  knots = stats::quantile(c(X), probs = c(seq(0, 1, 1/(num_knots - 1))))
  
  
  res_signal = list()
  
  for(signal_i in 1:5){
    res_iter= array(0, c(1,6))
    colnames(res_iter) = c("cons_diff","fun_diff","ise","b^2","||sum_mj||^2","||m||^2")
    for(iter in idx){
      BB = result[[iter]][[signal_i]]
      b0 = result[[iter]][[signal_i+5]]
      
      #the entry is \Vert demean(m_j) - demean(\hat{m}_j) \Vert^2
      fun_diff_ten = functionpart_diff_ten(BB, case=signal_i, precision=1000, knots = knots, order = 4)
      #the entry is \int \hat{m}_j  
      fhat_intten = fhat_ten_int(BB, case=signal_i, precision = 1000,knots=knots, order = 4)
      #the entry is \int m_j (dose not include the intercept)
      ftrue_intten = ftrue_mean_ten(case=signal_i)
      # sum(fhat_intten)+b0 is \hat{intercept}, 1 - sum(ftrue_intten) is true intercept
      cons_diff = (sum(fhat_intten)+b0 - 1 - sum(ftrue_intten))^2
      # \sum_j \Vert demean(m_j) - demean(\hat{m}_j) \Vert^2
      fun_diff = sum(fun_diff_ten)
      # the entry is \int m_j^2 (dose not include the intercept) 
      ftrue_square_intten = ftrue_demean_square_int_ten(case=signal_i)
      # sum(ftrue_square_intten) is \sum \int m_j^2, sum(ftrue_intten) + 1 is the intercept
      ftrue_norm_squre = sum(ftrue_square_intten) + (sum(ftrue_intten) + 1)^2
      
      true_b_squre = (sum(ftrue_intten) + 1)^2
      true_sum_mj = sum(ftrue_square_intten)
      ise_res = cons_diff + fun_diff
      res_iter[1,] = c(cons_diff, fun_diff, ise_res, true_b_squre, true_sum_mj, ftrue_norm_squre)
    }
    res_signal[[signal_i]] =  res_iter
  }
  
  n_seq_list[[samplesize_i]] = res_signal
  
  }
  return(n_seq_list)
}

sfSource('/data/yzhou/BNTR/FullSimPaper/CollectISEs/additional_utility_parallel.R')

sfClusterSetupRNG()

result_tem <- sfLapply(1:50, wrapper)


result = list()
n_name_seq = c("n=250", "n=500", "n=750", "n=1000")
signal_name_seq= c("case1","case2","case3","case4","case5")
for(n_name_i in 1:4){
n_name = n_name_seq[n_name_i] 
result_signal = list()
for(signal_i in 1:5){
  signal_name = signal_name_seq[signal_i]
  ise_meta = matrix(0, 50, 6)
  for(iter in 1:50){
    ise_meta[iter,] = result_tem[[iter]][[n_name_i]][[signal_i]]
  }
  result_signal[[signal_name]] = ise_meta
}
result[[n_name]] = result_signal
}



sfStop()
filename <- str_c("ise_meta_tnorm_n250_500_750_1000.Rdata")
save(result,file = filename)

















