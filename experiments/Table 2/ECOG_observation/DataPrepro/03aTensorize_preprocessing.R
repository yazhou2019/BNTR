#after the filter and wavenet transformation step, processe 

chall_time_frequen=read.csv('chall_time_frenquen_1.csv',header = F)

#y=read.table('y.csv',col.names = F)

X_input=array(0,c(64,10,10,10000))

for(j in 1:64){
  print(j)
  for(i in 1:10000){
    print(i)
    
    Matrix_new=matrix(as.numeric(chall_time_frequen[(i-1)*64+j,1:100]),10,10)
    X_input[j,,,i]=Matrix_new
    
    
  }
}

set.seed(2018)
id=sample(10000,10000)
idtrain=id[1:4000]

idtest=id[4001:10000]


X_train=array(0,c(64,10,10,4000))
X_test=array(0,c(64,10,10,6000))

for(i in 1:4000){
  X_train[,,,i]=X_input[,,,idtrain[i]]
}
for(i in 1:6000){
  X_test[,,,i]=X_input[,,,idtest[i]]
}

save(X_train,file="X_train_1.Rdata")
save(X_test,file="X_test_1.Rdata")

