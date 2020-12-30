# load the true coefficient tensor (matrix)
data("X_horse")
BB <- X_horse

set.seed(2019)

# sample size of the training and test set
n_train <- 400
n_test <- 100

# signal level
signal_level = 0.5


# the tuning parameters have been tuned by authors though validation; see the paper and test_validation.R
# rank used in the algorithm
rank <- 4
# the tuning parameters of the elastic net
lambda_1 <- 5
lambda_2 <- 0.5  # 0 corresponds to lasso; 1 corresponds to ridge
# use the initial point tuned from the validation (1) or not (0, sequential warmstart)
input_initial <- 0


# generate the training data
X_train = array(runif(prod(c(dim(BB), n_train)), 0, 1), c(dim(BB), n_train))
# broadcated procedure
BroX_train = X_train + 0.6 * sin(2 * pi * (X_train - 0.5)^2)
y_train = 1 + crossprod(matrix(BroX_train, c(prod(dim(BB)), n_train)), as.vector(BB)) + signal_level * rnorm(n_train)


# generate the test data
X_test = array(runif(prod(c(dim(BB), n_test)), 0, 1), c(dim(BB), n_test))
# broadcated procedure
BroX_test = X_test + 0.6 * sin(2 * pi * (X_test - 0.5)^2)
y_test = 1 + crossprod(matrix(BroX_test, c(prod(dim(BB)), n_test)), as.vector(BB)) + signal_level * rnorm(n_test)


# Transform to truncated power basis
num_knots <- 5
knots = quantile(c(X_train), probs = c(seq(0, 1, 1/(num_knots - 1))))
tildePhiX_train = tildePhiX_trans(X_train, knots)



# BroadcasTR
if (input_initial == 1) {
  data("initial_point")
  beta0 <- initial_point$beta0
  B0 <- initial_point$B0
  warmstart <- 0
} else {
  beta0 <- NA
  B0 <- NA
  warmstart <- 1
}
res = broadcasted_sparsetenreg(tildePhiX_train, y_train, r = rank, lambda = lambda_1, alpha = lambda_2, warmstart = warmstart, beta0 = beta0, B0=B0, Replicates=1)

# norm tensor
normtensor <- fhatnorm_ten(full_R(res$beta), knots)

# plot the true and estimated region
par(mar = c(0, 0, 0, 0), mfcol = c(1, 2), mai = c(0.2, 0.4, 0.4, 0.2))
plot(c(1, dim(BB)[1]), c(1, dim(BB)[2]), xaxt = "n", yaxt = "n", type = "n")
rasterImage(BB, 0, 0, dim(BB)[1] + 1, dim(BB)[2] + 1)
mtext("True coefficient tensor (matrix)", 3, line = 0.2)
plot(c(1, dim(BB)[1]), c(1, dim(BB)[2]), xaxt = "n", yaxt = "n", type = "n")
mtext("Estimated norm tensor (matrix)", 3, line = 0.2)
rasterImage(normtensor, 0, 0, dim(BB)[1] + 1, dim(BB)[2] + 1)


# prediction on the test set
tildePhiX_test = tildePhiX_trans(X_test, knots)
y_pre = res$beta0 + crossprod(matrix(tildePhiX_test, c(prod(dim(full_R(res$beta))), n_test)), as.vector(full_R(res$beta)))



# prediction performance
cat("The prediction performance in these tunning parameters: \n", "MSPE =", sum((y_test - y_pre)^2)/n_test, "\n")
