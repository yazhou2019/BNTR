library(divest)
library(RNifti)
library(BNTR)
library(rTensor)
library(Matrix)
library(glmnet)
p=c(80,80,3)


current = "~/Desktop/ADNI/PET/ADNI2"
setwd(current)
namelist=dir()
X_AD = array(0,c(length(namelist),prod(p)))
y_AD_name = c()
for(i in 1:length(namelist)){
  
  path = stringr::str_c(current,"/",namelist[i])
  
  restem=try(readDicom(path = path, subset = NULL, flipY = TRUE, crop = FALSE,
                       forceStack = FALSE, verbosity = 0L, labelFormat = "T%t_N%n_S%s",
                       depth = 5L, interactive = FALSE))
  if('try-error' %in% class(restem)){
    next
  }else{
    res=restem
  }
  
  ooo = res[[length(res)]]
  ooo=as.character(ooo)
  y_AD_name[i] = ooo
  ooo=strsplit(ooo,split="_")[[1]][[4]]
  #ooo=strsplit(ooo,split="N")[[1]][[2]]
  #print(strsplit(ooo,split="_")[[1]][[2]])
  
  res=asNifti(res[[length(res)]])
  res=as.array(res)
  res=as.numeric(res)
  res=array(res,c(160,160,96))
  res=res - sum(res)/(160*160*96)
  
  res=res[,,51:53]
  res=resize_array(res,p)
  
  X_AD[i,]=c(res) - mean(c(res))
  #  y_AD_name[i] = ooo
}





current = "~/Desktop/ADNI/CNDATA/ADNI"
setwd(current)
namelist=dir()
X_CN = array(0,c(length(namelist),prod(p)))
y_CN_name = c()
for(i in 1:length(namelist)){
  
  path = stringr::str_c(current,"/",namelist[i])
  
  restem=try(readDicom(path = path, subset = NULL, flipY = TRUE, crop = FALSE,
                       forceStack = FALSE, verbosity = 0L, labelFormat = "T%t_N%n_S%s",
                       depth = 5L, interactive = FALSE))
  if('try-error' %in% class(restem)){
    next
  }else{
    res=restem
  }
  
  ooo = res[[length(res)]]
  ooo=as.character(ooo)
  y_CN_name[i] = ooo
  ooo=strsplit(ooo,split="_")[[1]][[4]]
  #ooo=strsplit(ooo,split="N")[[1]][[2]]
  
  res=asNifti(res[[length(res)]])
  res=as.array(res)
  res=as.numeric(res)
  res=array(res,c(160,160,96))
  res=res - sum(res)/(160*160*96)
  
  res=res[,,51:53]
  res=resize_array(res,p)
  X_CN[i,]=c(res) - mean(c(res))
  #  y_CN_name[i]=ooo
}

X = rbind(cbind(X_AD),cbind(X_CN))
y = c(rep(1,dim(X_AD)[1]), rep(-1,dim(X_CN)[1]))
y_names = c(y_AD_name, y_CN_name)

#datawithname = list(X_all = X, y_all = y, names = y_names)
#save(file = "datawithname80803.Rdata", datawithname)
#save(file = "imagenames.Rdata", y_names)
#只有943个，第944个是0

for(j in 1:943){
  if(sum(X[j,]==rep(0,prod(p)))==prod(p)){print(j)}
}