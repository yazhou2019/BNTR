library(R.matlab)


load("./data/X_train_1.Rdata")
load("./data/X_test_1.Rdata")
load("./data/y_all.Rdata")


writeMat("X_train_1_MATLAB.mat", X_train=X_train)
writeMat("X_test_1_MATLAB.mat", X_test=X_test)
writeMat("y_all_MATLAB.mat", y_all=y_all)



