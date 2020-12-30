data("X_horse")
BB <- X_horse

set.seed(2019)

# sample size of the training and test set
n_train <- 400
n_test <- 100

# signal level
signal_level = 0.5


# training data
X_train = array(runif(prod(c(dim(BB), n_train)), 0, 1), c(dim(BB), n_train))
# broadcated procedure
BroX_train = X_train + 0.6 * sin(2 * pi * (X_train - 0.5)^2)
y_train = 1 + crossprod(matrix(BroX_train, c(prod(dim(BB)), n_train)), as.vector(BB)) + signal_level * rnorm(n_train)

# the test data
X_test = array(runif(prod(c(dim(BB), n_test)), 0, 1), c(dim(BB), n_test))
BroX_test = X_test + 0.6 * sin(2 * pi * (X_test - 0.5)^2)
y_test = 1 + crossprod(matrix(BroX_test, c(prod(dim(BB)), n_test)), as.vector(BB)) + signal_level * rnorm(n_test)

set.seed(1224) # used for validation method
n_vali = 0.25 * n_train
# generate the test data
X_vali = array(runif(prod(c(dim(BB), n_vali)), 0, 1), c(dim(BB), n_test))
BroX_vali = X_vali + 0.6 * sin(2 * pi * (X_vali - 0.5)^2)
y_vali = 1 + crossprod(matrix(BroX_vali, c(prod(dim(BB)), n_vali)), as.vector(BB)) + signal_level * rnorm(n_vali)


R=c(1,2,3,4,5)
alpha=c(0,0.5,1)
lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)

res_valiandtest <- validation_broadcasted_sparsetenreg(R,alpha,lambda,X_train,y_train,X_vali,y_vali,X_test,y_test, num_knots=5, order=4)


