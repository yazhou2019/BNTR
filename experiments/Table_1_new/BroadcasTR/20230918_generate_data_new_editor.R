#generate the simulation data

source('./ComponentsNonLin/functions_needed.R')
source('./ComponentsNonLin/sequential_warmstart.R')
source('./ComponentsNonLin/broadcasted_sparsetenreg.R')
source('./ComponentsNonLin/validation_broadcasted_sparsetenreg.R')

source("./SimDataGeneration_He/20230224_get_outputs.R")
source("./SimDataGeneration_He/20230224_get_inputs.R")
source("./SimDataGeneration_He/20230224_get_BB.R")
source("./SimDataGeneration_He/20230224_get_function.R")


n_all = 1000

set.seed(2022)

v = list()
for(i in 1:5){
  v[[i]]=1
}


BBprod = get_BB()
data_all = list()
for (iter in 1:50){
  set.seed(iter)
  X_data = get_inputs(n=n_all, p=c(64,64), distri = "uniform")
  y_all = get_outputs(v, BBprod, X_data)
  data_all[[iter]] = list(X=X_data, y=y_all)
}

for (iter in 1:50){
  print("itertation")
  print(iter)
  for (CASE in 1:5){
    data_all[[iter]]$y[[CASE]] = data_all[[iter]]$y[[CASE]] + rnorm(n_all, 0.1*var(data_all[[1]]$y[[CASE]])^0.5)
  }
}

data_all[[51]] = BBprod
data_all[[52]] = v


filename <- str_c("Simul_n", n_all, "_rep50_final_fix_new_editor.Rdata")
setwd("./SimResults")
save(data_all,file = filename)
setwd("../")


