#save the output

y=read.table('y_1.csv',col.names = F)
Y_big=y

save(Y_big,file="y_output_1.Rdata")


set.seed(2018)
id=sample(10000,10000)
idtrain=id[1:4000]
idtest=id[4001:10000]

y_test=Y_big[idtest]
y_train=Y_big[idtrain]

y_all=c()
y_all[1:4000]=y_train
y_all[4001:10000]=y_test

save(y_all,file="y_all.Rdata")

#save the output

#for big data set
#y=read.table('y_1_big.csv',col.names = F)
#Y_big=y

#save(Y_big,file="y_output_1_big.Rdata")
