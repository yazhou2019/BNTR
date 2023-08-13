library(readr)
library(stringi)
library(stringr)

# data(mmse)
# mmsemach(mmse$RID, mmse$USERDATE, mmse$MMSCORE)
mmsemach <- function(rid,dtime,mmse_value){
  n_all = length(dtime)
  dtime = as.character(dtime)
  dtime_rid = c()
  for(i in 1:n_all){
    dt = strsplit(dtime[i],split="-")[[1]]
    year = dt[1]
   #month =  as.numeric(dt[2])
   #day =  as.numeric(dt[3])
    
    # if(month <= 9){
    #   month = str_c("0",as.character(month))
    # }else{
    #   month = as.character(month)
    # }
    # if(day <= 9){
    #   day = str_c("0",as.character(day))
    # }else{
    #   day = as.character(day)
    # }
    
    #dt = str_c(year, month, day)
    dt = year
    dtime_rid[i] = str_c(dt, "_", as.character(rid[i]))
  }
 res = cbind(dtime_rid, mmse_value) 
 return(res)
}
  



namemach <- function(sb,dtime){
  dtime_sb = c()
  
  n_all = length(dtime)
  for(i in 1:n_all){
    dt = strsplit(dtime[i],split="/")[[1]]
    year = dt[3]
    month =  as.numeric(dt[1])
    day =  as.numeric(dt[2])
    
    if(month <= 9){
      month = str_c("0",as.character(month))
    }else{
      month = as.character(month)
    }
    if(day <= 9){
      day = str_c("0",as.character(day))
    }else{
      day = as.character(day)
    }
    
    dt = str_c(year, month, day)
    
    dtime_sb[i] = str_c(dt, "_",sb[i])
  }
  return(dtime_sb)
}

namesub <- function(y){
  n_all = length(y)
  for(i in 1:n_all){
    if(is.na(y[i])!=1){
      tem = strsplit(y[i],split="_")[[1]]
      dt = substr(tem[1],2,9)
      sbname = str_c('_',substr(tem[2],2,4),'_',tem[3],'_',tem[4])
      y[i] = str_c(dt, sbname)
    }
  }
  return(y)
}


sexage_filter <- function(sexage){
  sex = sexage[1]
  age =sexage[2]
  descirp = sexage[3]
  mmse_value = as.numeric(sexage[6])
  if(sex=="M"){
    sex = 1
  }else{
    sex = 2
  }
  age = as.numeric(age)
  return(c(sex,age,as.numeric(descirp),mmse_value))
}

description_filter = function(descrip){
  n_all = length(descrip)
  isornot = c()
  for(i in 1:n_all){
    #if(descrip[i] == "Coreg, Avg, Standardized Image and Voxel Size"){
    if(descrip[i] == "AV45 Coreg, Avg, Standardized Image and Voxel Size"){
      isornot[i] = 1
    }else{
      isornot[i] = 0
    }
  }
  return(isornot)
}

AD_12_25_2020 <- read_csv("./ADNI/AD_12_25_2020.csv")
AD_file = namemach(sb = AD_12_25_2020$Subject, dtime = AD_12_25_2020$`Acq Date`)
AD_file = cbind(AD_file, AD_12_25_2020$Sex, AD_12_25_2020$Age,description_filter(AD_12_25_2020$Description))
CN_12_25_2020 <- read_csv("./ADNI/CN_12_25_2020.csv")
CN_file = namemach(sb = CN_12_25_2020$Subject, dtime = CN_12_25_2020$`Acq Date`)
CN_file = cbind(CN_file, CN_12_25_2020$Sex, CN_12_25_2020$Age,description_filter(CN_12_25_2020$Description))
ADCN_file = rbind(AD_file, CN_file)

temname = c()
for(i in 1:dim(ADCN_file)[1]){
 tem = strsplit(ADCN_file[i,1],split="_")[[1]]
 tem4 = tem[4]
 tem4 = as.numeric(tem4)
 tem4 = as.character(tem4)
 temname[i] = str_c(substr(tem[1],1,4), "_", tem4) 
}
ADCN_file = cbind(ADCN_file, temname)


mmse_file = mmsemach(mmse$RID, mmse$USERDATE, mmse$MMSCORE)


tem_mmse = cbind(rep(NA, dim(ADCN_file)[1]),rep(NA, dim(ADCN_file)[1]))
for(i in 1:dim(ADCN_file)[1]){
  tem = ADCN_file[i,5]
  for(j in 1:dim(mmse_file)[1]){
    if(tem == mmse_file[j,1]){
      tem_mmse[i,]=mmse_file[j,]
    }
  }
}
ADCN_file = cbind(ADCN_file, tem_mmse)



load("./datawithname40403.Rdata")

X=datawithname$X_all
y=datawithname$y_all
y_names = datawithname$names
n_all = length(y_names)
new_y_names = namesub(y_names) 
sex_age = matrix(0, n_all, 4)
for(i in 1:n_all){
    index = max(which(ADCN_file[,1]==new_y_names[i]))
    sex_age[i,]=sexage_filter(ADCN_file[index,2:7])
}

sex_age_des_mmse=sex_age

data = list(X_all = X, y_all = y, names= y_names, sex_age=sex_age)
save(file="data40403.Rdata", data_all)
