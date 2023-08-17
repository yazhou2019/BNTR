# there are some functions that are needed for the algorithm

# generate \phi(X)
tildePhiX_trans <- function(X_sample, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4) {

    num_knots <- length(knots)

    # number of basis
    K <- num_knots - 2 + order - 1
    if (K == 1)
        return(X_sample)

    p1 <- dim(X_sample)
    n <- p1[length(p1)]
    p <- p1[-length(p1)]
    BX <- array(0, c(p, K, n))
    str_slice <- c()
    for (j in 1:(length(p1) + 1)) {
        if (j == (length(p) + 1)) {
            str_slice <- c(str_slice, "k")
        } else {
            str_slice <- c(str_slice, ",")
        }
    }
    str_slice <- paste(str_slice, collapse = " ")
    str_text1 <- paste("BX[", str_slice, "]=X_sample^k", collapse = " ")
    str_text2 <- paste("BX[", str_slice, "]=((X_sample - knots[k-order+1])*((X_sample - knots[k-order+1])>0))^(order - 1)",
        collapse = " ")

    for (k in 1:K) {
        if (k < order) {
            eval(parse(text = str_text1))
        } else {
            eval(parse(text = str_text2))
        }
    }
    return(BX)
}











linearpart_penalty_function = function(y, X, beta, beta0, lambda) {
    X = X@data
    D = full_R(beta)
    p1 = dim(X)
    n = p1[length(p1)]
    linearpart = c()
    for (i in 1:n) {
        linearpart[i] = sum(c(D) * c(as.array(X[, , , i]))) + beta0[1]
    }

    l = loglike(y, linearpart, family = "gaussian")
    P = l + n * lambda * penaltyvalue(beta, alpha = 0, manifold = 0)
    L = l
    test = list(L = L, lossP = P)
    return(test)
}


loglike <- function(y, linearpart, family = "binomial") {
    if (family == "binomial") {
        dev = crossprod(y, linearpart) - sum(log(1 + exp(linearpart[linearpart <= 600]))) - sum(linearpart[linearpart >
            600])
        return(dev)
    }
    # guassian L2 loss
    if (family == "gaussian") {
        return(sum((y - linearpart)^2))
    }
}


# transform the CP components to a tensor
full_R <- function(Beta) {
    d = length(Beta)
    p = matrix(0, d, 1)
    r = dim(Beta[[1]])[2]
    for (i in 1:d) {
        p[i] = dim(Beta[[i]])[1]
    }
    B = array(0, c(p))
    for (j in 1:r) {
        B1 = Beta[[1]][, j]
        for (i in 2:d) {
            B1 = B1 %o% Beta[[i]][, j]
        }
        B = B + B1
    }
    return(B)
}

arrange <- function(beta, a, manifoldfix = 1) {
    # the old version r and d are inter-changed. r here is D+1 r=length(beta)
    if (manifoldfix == 1) {
        d = length(beta) - 1
    } else {
        d = length(beta)
    }
    # d here is r d=dim(beta[[1]])[2]
    r = dim(beta[[1]])[2]
    if (a >= 0) {
        a = a^(1/d)
        for (i in 1:d) {
            beta[[i]] = beta[[i]] * a
        }
    } else {
        a = (-a)^(1/d)
        for (i in 1:d) {
            if (i == 1) {
                beta[[i]] = -beta[[i]] * a
            } else {
                beta[[i]] = beta[[i]] * a
            }
        }

    }

    return(beta)
}

# the rescale factor for elastic this is used for Newton method in rescale strategy#### f_d^\prime(x)
elastic_gradient_d = function(betad, alpha, lambda) {
    betadl1 = sum(abs(betad))
    betadl2square = sum((betad)^2)

    numerator = 0.5 * 4 * (1 - alpha) * betadl2square
    denominator1 = -alpha * betadl1 + sqrt(alpha^2 * betadl1^2 + 4 * lambda * (1 - alpha) * betadl2square)
    # print(denominator1) print('this') print(betadl1) print(lambda) print(betadl2square) print('that')
    denominator2 = sqrt(alpha^2 * betadl1^2 + 4 * lambda * (1 - alpha) * betadl2square)
    denominator = denominator1 * denominator2

    gradient_for_d = numerator/denominator
    return(gradient_for_d)

}

# f_d(x)
elastic_f_lambda_d = function(betad, alpha, lambda) {
    betadl1 = sum(abs(betad))
    betadl2square = sum((betad)^2)
    term1 = log(-alpha * betadl1 + sqrt(alpha^2 * betadl1^2 + 4 * lambda * (1 - alpha) * betadl2square))
    term2 = log(betadl2square)
    termd = term1 - term2
    return(termd)
}

# f/f^\prime
elastic_foverfprime = function(beta, alpha, lambda, d, rr) {
    numerator_f = 0
    denominator_f_prime = 0
    for (i in 1:d) {
        betad = beta[[i]][, rr]
        numerator_f = numerator_f + elastic_f_lambda_d(betad, alpha, lambda)
        denominator_f_prime = denominator_f_prime + elastic_gradient_d(betad, alpha, lambda)
    }
    numerator_f = numerator_f - d * log(2 - 2 * alpha)
    foverfprime = numerator_f/denominator_f_prime
    return(c(foverfprime, numerator_f))
}

# Newton's method for r return the lagrange factor
elastic_newton_r = function(beta, alpha, d, rr) {
    lambda = 1
    for (i in 1:d) {
        betad = beta[[i]][, rr]
        lambda = lambda * sum(betad^2)
    }
    lambda = lambda^(1/d)
    # the total iteration for newton is 10
    for (iter in 1:10) {
        lambda = lambda - elastic_foverfprime(beta, alpha, lambda, d, rr)[1]
        # print('hehe') print(elastic_foverfprime(beta,alpha,lambda,d,rr)[2]) later on add a if condition
    }
    return(lambda)
}


# here lambda is the lagrange, not the tunning parameter in the model return the rescale factor by Newton,
# which need to be adjusted.
elastic_rescale_factor_rough = function(betad, alpha, lambda) {
    betadl1 = sum(abs(betad))
    betadl2square = sum((betad)^2)
    # the rescale factor
    numerator = -alpha * betadl1 + sqrt(alpha^2 * betadl1^2 + 4 * lambda * (1 - alpha) * betadl2square)
    denominator = 2 * (1 - alpha) * betadl2square
    epsilond = numerator/denominator
    return(epsilond)
}


# return the rescale factor
elastic_rescale_factor = function(beta, alpha, d, rr) {
    epsilon_all = 1
    lambda = elastic_newton_r(beta, alpha, d, rr)
    # print('test') print(elastic_foverfprime(beta,alpha,lambda,d,rr)[2]) print('test')

    epsilon_d = c()
    for (i in 1:d) {
        betad = beta[[i]][, rr]
        epsilon_d[i] = elastic_rescale_factor_rough(betad, alpha, lambda)
        epsilon_all = epsilon_all * epsilon_d[i]
    }
    # print(epsilon_all) print(elastic_foverfprime(beta,alpha,lambda,d,rr)[2])
    tuningto1 = (1/epsilon_all)^(1/d)
    # print(tuningto1) print(epsilon_d)
    for (i in 1:d) {
        epsilon_d[i] = epsilon_d[i] * tuningto1
    }
    # print(epsilon_d)
    return(epsilon_d)
}



# if beyond the limitation, then, no scale all=1 means rescale all components, otherwise, rescale the former
# d-1 components. l2=1 means
rescale <- function(beta, rescale_all = 1, rescale_l2 = 1, elastic = 0, alpha = 0) {
    if (rescale_all == 1) {
        d = length(beta)
    } else {
        # only rescale CP coefficients
        d = length(beta) - 1
    }
    r = dim(beta[[1]])[2]
    lambda_r = rep(1, r)
    # in some step, part of beta has been chanted elastic rescale strategy is also included in the model
    beta_org = beta

    for (i in 1:r) {
        for (ddd in 1:d) {
            if (rescale_l2 == 1) {
                lambda_ir = sum(beta[[ddd]][, i]^2)^(0.5)
            } else {
                lambda_ir = sum(abs(beta[[ddd]][, i]))
            }
            # to avoid numerical problem
            if (is.nan(lambda_ir) == 1) {
                return(beta)
            }
            # to avoid numerical problem
            if (lambda_ir == 0) {
                return(beta_org)
            }
            lambda_r[i] = lambda_r[i] * lambda_ir
            beta[[ddd]][, i] = beta[[ddd]][, i]/lambda_ir
        }
    }
    lambda_r = lambda_r^{
        1/d
    }
    for (i in 1:r) {
        for (ddd in 1:d) {
            beta[[ddd]][, i] = beta[[ddd]][, i] * lambda_r[i]
        }
    }


    beta_elsatic = beta
    if (elastic != 0 && alpha != 1) {
        for (rr in 1:r) {
            rescaleforr = elastic_rescale_factor(beta, alpha, d, rr)
            for (ddd in 1:d) {
                beta_elsatic[[ddd]][, rr] = beta[[ddd]][, rr] * rescaleforr[ddd]
            }

        }

        # to avoid numeraical problems
        if (is.na(sum(rescaleforr)) != 1) {
            # print('OK') print(rescaleforr)
            beta = beta_elsatic
        }
    }





    return(beta)
}



predict_BroadcasTR <- function(b0, BB, X, y=NA, family = "gaussian", returnyhat = 0){
  D1=dim(X)
  n=D1[length(D1)]
  if(family=="gaussian"){
  y_hat = b0 + crossprod(matrix(X, c(prod(dim(BB)), n)), as.vector(full_R(BB)))
  if(returnyhat!=0){
    return(y_hat)
  }else{
    MPSE = sum((y_hat - y)^2)/n
    return(MPSE)
  }
  }

}

prediction_function_nonlinear <- function(b0, BB, X, y, family = "gaussian", yyyreturn = 0) {
    d = length(dim(X)) - 1
    if (family == "gaussian") {
        n = length(y)
        BB = rTensor::as.tensor(BB)
        if (class(X)[1] == "Tensor") {
        } else {
            X = rTensor::as.tensor(X)
        }
        yhat = c()



        if (d == 4) {
            BBB = array(0, c(1, dim(BB)))
            BBB[1, , , , ] = BB@data
            yhat = t(ctprod(BBB, X@data, 4)) + b0
        }

        if (d == 3) {
            BBB = array(0, c(1, dim(BB)))
            BBB[1, , , ] = BB@data
            yhat = t(ctprod(BBB, X@data, 3)) + b0
        }




        rm(X)
        gc()
        MSE = norm(as.matrix(y) - yhat, "f")^2/n
        if (yyyreturn == 0) {
            return(MSE)
        } else {
            test = list(MSE = MSE, yhat = yhat)
        }
    }
}
prediction_function_linear <- function(b0, BB, X, y, family = "gaussian") {
    if (family == "gaussian") {
        n = length(y)
        if (class(BB)[1] == "Tensor") {
        } else {
            BB = rTensor::as.tensor(BB)
        }


        if (class(X)[1] == "Tensor") {
        } else {
            X = rTensor::as.tensor(X)
        }

        yhat = c()
        d = length(dim(X)) - 1

        if (d == 3) {
            for (i in 1:n) {
                yhat[i] = rTensor::innerProd(BB, X[, , , i]) + b0
            }
        }
        if (d == 2) {
            for (i in 1:n) {
                yhat[i] = rTensor::innerProd(BB, X[, , i]) + b0
            }
        }
        # RPE
        MSE = norm(as.matrix(y) - yhat, "f")^2/n
        return(MSE)
    }
}

rlm <- function(X, y, uk) {
    return(0)
}
# X_sample=array(rnorm(10*5*6*20),c(10,5,6,20))


# optimization on the oblique manifold

grad <- function(Y, Xk_vec, D, K, r) {

    nablaLvec = (-2 * t(D) %*% Y + 2 * t(D) %*% D %*% Xk_vec)
    nablaL = matrix(nablaLvec, K, r)
    Xk = matrix(Xk_vec, K, r)

    # print(dim(Xk)) print(dim(nablaL))
    if (r == 1) {
        gradf = nablaL - Xk %*% (t(Xk) %*% nablaL)
    } else {
        gradf = nablaL - Xk %*% diag(diag(t(Xk) %*% nablaL))
    }
    return(gradf)
}
Retraction <- function(grad, alpha, Xk) {
    Xk1 = Xk + alpha * (-grad)
    r = dim(Xk1)[2]
    if (r == 1) {
        Xk1 = Xk1 %*% (t(Xk1) %*% Xk1)^(-0.5)
    } else {
        Xk1 = Xk1 %*% diag((diag(t(Xk1) %*% Xk1))^(-0.5))
    }
    return(Xk1)
}

falphaL = function(D, X, Y) {
    return(sum((Y - D %*% c(X))^2))
}

update_alpha = function(fast = 0, alphafast = 0, Y, BD1, D, K, r, maxupdate_out = 1, maxupdate_in = 40, scalefactor = 0.5,
    thresholdformani = 1e-10, guess = 0) {
    X0 = BD1
    # print(dim(D)) print(dim(X0)) print(dim(Y))
    f0 = falphaL(D, X0, Y)
    Xl = X0
    fl = f0
    fcontrol = f0

    gradX0 = grad(Y, c(X0), D, K, r)
    gradXl = gradX0

    # the outer loop, both for fast or not fast
    for (l in 1:maxupdate_out) {
        # print(l)
        norm2gradXl = norm(gradXl, "f")^2
        # print(norm2gradXl)

        # the basic method to choose alpha, guess or 0.5^l*alpha
        if (guess == 1) {
            alpha = 1/(norm2gradXl^0.5)
            print(alpha)
            print(fl)
        } else {
            alpha = 2
        }


        # this will be used for those cases that do not satisfy Armijo
        Xmay = Xl
        fmay = fl
        # this is used for those cases that the monifold dose not move fcontrol=fl

        # index for Armijo criterion
        ArmijoControl = 0

        # The fast setting
        ArmijoControlfast = 0
        alphafastuse = alphafast

        # this two indexes are used to
        controlindexmore = 0
        controlindexless = 0

        # The fast setting
        if (fast == 1) {
            # this loop is to find a small step size
            for (i in 1:maxupdate_in) {
                Xl1 = Retraction(grad = gradXl, alpha = alphafastuse, Xk = Xl)
                fl1 = falphaL(D = D, X = Xl1, Y = Y)

                if ((fl1 - fl) <= -0.5 * alphafastuse * norm2gradXl) {
                  if (controlindexless == 1) {
                    # print('Armijofast')
                    ArmijoControlfast = 1
                    Xl = Xl1
                    fl = fl1
                    gradXl = grad(Y, c(Xl), D, K, r)
                    break
                  }
                  controlindexmore = 1
                  alphafastuse = 2 * alphafastuse

                } else {
                  controlindexless = 1
                  alphafastuse = 0.5 * alphafastuse

                }
            }
        }

        if (ArmijoControlfast == 1) {
            alphafast = alphafastuse
            break
        }


        # the old and basic one
        for (i in 1:maxupdate_in) {
            if (guess == 1) {
                # this guess may need to be adjusted
                alpha = alpha * scalefactor
            } else {
                alpha = alpha * scalefactor
            }


            # print(alpha)
            Xl1 = Retraction(grad = gradXl, alpha = alpha, Xk = Xl)
            # print(Xl1) print(dim(D)) print(dim(Xl1))
            fl1 = falphaL(D = D, X = Xl1, Y = Y)
            # print(fl1)
            if ((fl1 - fl) <= (-0.5 * alpha * norm2gradXl)) {
                ArmijoControl = 1
                Xl = Xl1
                fl = fl1
                # print(fl1)
                alphafast = alpha
                gradXl = grad(Y, c(Xl), D, K, r)
                break
            }

            # if decline, we can accept it
            if (fl1 < fl) {
                Xmay = Xl1
                fmay = fl1
                alphafast = alpha
            }

        }


        # print(fl1-fl) print(-0.5*alpha*norm2gradXl)


        # not use Armijo step-size rule
        if (ArmijoControl == 1) {
            # print('Armijo')
        } else {
            if (fmay >= fcontrol) {

                # print('manifold optimization cannot move')
                Notmove = list(Xl = X0, alphafast = 1)
                return(Notmove)
                break
            } else {
                # print('Not Armijo') print(fl) print(fmay)
                Xl = Xmay
                fl = fmay
                gradXl = grad(Y, c(Xl), D, K, r)
            }
            #
        }
        #

        if (max(abs(gradXl)) <= thresholdformani * (1 + max(abs(gradX0)))) {
            print("good")
            break
        }

    }

    if (fl > f0) {
        result = list(Xl = X0, alphafast = alphafast, fl = f0)
    } else {
        result = list(Xl = Xl, alphafast = alphafast, fl = fl)
    }
    return(result)
}


# only used when elastic nett penalty, later may works for other penalty.
penaltyvalue = function(beta, alpha, manifold = 1, QCQP = 0) {
    if (manifold == 1 || QCQP == 1) {
        D = length(beta) - 1
    } else {
        D = length(beta)
    }


    value = 0
    value1 = 0
    value2 = 0
    if (alpha == 0) {
        for (i in 1:D) {
            value = value + sum((beta[[i]])^2)
        }
    } else {
        for (i in 1:D) {
            value1 = value1 + sum(abs(beta[[i]]))
        }
        value1 = 2 * value1
        for (i in 1:D) {
            value2 = value2 + sum((beta[[i]])^2)
        }
        value2 = value2

        value = value1 + value2

    }
    return(value)
}

# for debug
distancefunction = function(linearpartlist) {
    d = length(linearpartlist)
    distancematrix = matrix(0, d, d)
    for (i in 1:d) {
        for (j in i:d) {
            distancematrix[i, j] = (sum((linearpartlist[[i]] - linearpartlist[[j]])^2))^0.5
        }
    }
    return(distancematrix)
}


# follow MATLAB TensorReg
linspace = function(d1, d2, n) {
    n = floor(n)
    n1 = n - 1
    c = (d2 - d1) * (n1 - 1)
    if (c == Inf) {
        if ((d2 - d1) == Inf) {
            y = d1 + (d2/n1) * c(0:n1) - d1/n1 * (0:n1)
        } else {
            y = d1 + c(0:n1) * ((d2 - d1)/n1)
        }
    } else {
        y = d1 + c(0:n1) * (d2 - d1)/n1
    }
    if (is.na(y[1])) {
    } else {
        if (d1 == d2) {
            y = d1
        } else {
            y[1] = d1
            y[length(y)] = d2
        }
    }
    return(y)
}

# follow MATLAB TensorReg
resize_array = function(Xarray, pobject) {
    p = dim(Xarray)
    D = length(p)
    d = length(pobject)
    U = list()
    Xarray_new = Xarray
    if (class(Xarray_new)[1] == "Tensor") {
        tensor_indicator = 1
    } else {
        tensor_indicator = 0
        Xarray_new = rTensor::as.tensor(Xarray_new)
    }
    for (d in 1:D) {
        xi = linspace(1, p[d], pobject[d])
        j1 = floor(xi)
        j2 = j1 + 1
        w1 = j2 - xi
        w2 = 1 - w1
        j2[length(j2)] = p[d]
        Ud = as.matrix(Matrix::spMatrix(pobject[d], p[d], c(1:pobject[d], 1:pobject[d]), c(j1, j2), c(w1, w2)))
        Xarray_new = rTensor::ttm(Xarray_new, Ud, m = d)
    }
    # return the same class as Xarray
    if (tensor_indicator == 1) {
        return(Xarray_new)
    } else {
        return(Xarray_new@data)
    }
}





# dual ascent method for alpha
dualascentforalpha <- function(y, DX, alphaold, C = 1, Maxiterfordual = 1000, TolFundualascent = 1e-07, MaxIterInner = 40) {
    r = dim(alphaold)[2]
    K = dim(alphaold)[1]
    initialvalue = c()
    h_alpha_intial = c()
    initialsacale = 2
    # initial value

    for (i in 1:10) {
        initialpointtest = MASS::ginv(t(DX) %*% DX + (initialsacale^i) * diag(r * K)) %*% t(DX) %*% y

        for (j in 1:r) {
            h_alpha_intial[j] = sum(initialpointtest[((j - 1) * K + 1):(j * K)]^2) - C
        }
        initialvalue[i] = abs(sum(initialpointtest * initialpointtest) - r * C)
        initialindex = sum(h_alpha_intial <= 0)
        if (initialindex == r) {

        } else {
            initialvalue[i] = Inf
        }

        if (initialvalue[i] < 0.1) {
            break
        }
    }
    initial_index = which.min(initialvalue)
    # print(initial_index)

    lambda0 = rep(initialsacale^initial_index, r)
    lambda = lambda0

    n = dim(DX)[1]

    lambda_diag_vec = rep(1, r * K)
    for (i in 1:r) {
        lambda_diag_vec[((i - 1) * K + 1):(i * K)] = lambda[i] * lambda_diag_vec[((i - 1) * K + 1):(i * K)]
    }
    lambda_diag = diag(lambda_diag_vec)
    h_alpha = rep(0, r)


    test_lambda = list()
    test_L_x = c()
    test_f = c()
    test_alpha = list()

    test_halpha = list()

    for (iter in 1:Maxiterfordual) {

        if (iter == 1) {
            alphavec = MASS::ginv(t(DX) %*% DX + lambda_diag) %*% t(DX) %*% y
            for (i in 1:r) {
                h_alpha[i] = sum(alphavec[((i - 1) * K + 1):(i * K)]^2) - C
            }

            f_0 = sum((y - DX %*% alphavec)^2)
            L_x_lambda_0 = f_0 + sum(h_alpha * lambda)

            test_lambda[[1]] = lambda
            test_alpha[[1]] = alphavec
            test_L_x[1] = L_x_lambda_0

            test_f[1] = f_0

            test_halpha[[1]] = h_alpha

        } else {
            tk = 2
            L_x_lambda_tem_final = L_x_lambda_0
            indexarmijo = 0
            h_alpha_final = h_alpha
            lambda_final = lambda
            f0_final = f_0
            alphavec_final = alphavec

            for (j in 1:MaxIterInner) {
                tk = tk * 0.5

                for (i in 1:r) {
                  # h_alpha[i]=sum(alphavec[((i-1)*K+1):(i*K)]^2)-C
                  lambda[i] = max((lambda[i] + tk * h_alpha[i]), 0)
                }

                lambda_diag_vec = rep(1, r * K)
                for (i in 1:r) {
                  lambda_diag_vec[((i - 1) * K + 1):(i * K)] = lambda[i] * lambda_diag_vec[((i - 1) * K + 1):(i *
                    K)]
                }

                lambda_diag = diag(lambda_diag_vec)

                alphavec = MASS::ginv(t(DX) %*% DX + lambda_diag) %*% t(DX) %*% y

                f_0 = sum((y - DX %*% alphavec)^2)

                for (i in 1:r) {
                  h_alpha[i] = sum(alphavec[((i - 1) * K + 1):(i * K)]^2) - C
                  # lambda[i]=max((lambda[i]+tk*h_alpha[i]),0)
                }

                L_x_lambda_tem_before = f_0 + sum(h_alpha * lambda)
                if (L_x_lambda_tem_before - L_x_lambda_0 > (tk * 0.5 * sum(h_alpha^2))) {
                  # feasbile adjustment
                  if (sum(h_alpha <= 0) == r) {
                    indexarmijo = 1
                    L_x_lambda_tem = L_x_lambda_tem_before

                    break
                  } else {
                  }
                }

                if (L_x_lambda_tem_before > L_x_lambda_tem_final) {
                  L_x_lambda_tem_final = L_x_lambda_tem_before
                  h_alpha_final = h_alpha
                  lambda_final = lambda
                  f0_final = f_0
                  alphavec_final = alphavec
                }

            }

            if (indexarmijo == 1) {
                # print('hehe')
            } else {
                h_alpha = h_alpha_final
                lambda = lambda_final
                L_x_lambda_tem = L_x_lambda_tem_final
                f_0 = f0_final
                alphavec = alphavec_final
            }


            L_x_lambda_diff = L_x_lambda_tem - L_x_lambda_0
            L_x_lambda_0 = L_x_lambda_tem




            test_lambda[[iter]] = lambda
            test_alpha[[iter]] = alphavec


            test_halpha[[iter]] = h_alpha

            feasibleindex = sum(h_alpha <= (0))
            if (feasibleindex == r) {
                test_L_x[iter] = L_x_lambda_0
                test_f[iter] = f_0
            } else {
                test_L_x[iter] = -Inf
                test_f[iter] = Inf
            }

            if (max(lambda) == 0) {
                break
            }

            if (abs(L_x_lambda_diff) <= TolFundualascent * (L_x_lambda_0)) {
                if (iter >= 10) {
                  break
                } else {
                }
            }

        }
    }



    indexf = which.max(test_L_x)

    # indexf=which.min(test_f)
    index = indexf
    if (sum(test_halpha[[indexf]] <= 0) == r) {

        alphavec_final = test_alpha[[indexf]]
    } else {

        alphavec_final = c(alphaold)
        # alphavec_final=test_alpha[[indexf]]
    }

    res = list(test_alpha = test_alpha, alphavec_final = alphavec_final, test_L_x = test_L_x, test_f = test_f,
        test_lambda = test_lambda, test_halpha = test_halpha, index = index)

    return(res)
}






warmstartdowsize = function(warmstart = 3, r = 2, X_sample = NA, shrink_factor_scale = 5, B0 = NA, B0_list = NA,
    BB0 = NA, BB0_list = NA) {
    # From a lot of experiment, warmstart=3 is the best
    p1 = dim(X_sample)
    n = p1[length(p1)]
    p = p1[-length(p1)]
    K = p[length(p)]
    d = length(p)


    B0_list_reduce = B0_list
    B0_reduce = B0
    BB0_reduce = BB0
    BB0_list_reduce = BB0_list

    X_sample_reduce = X_sample

    shrink_factor = (n/shrink_factor_scale)/(r * sum(p))

    if (warmstart == 3) {
        # print(shrink_factor)
        if (shrink_factor < 1) {
            # print(shrink_factor) print(p) print(c(p[1:(d-1)]*shrink_factor,K)) heheda=c(p[1:(d-1)]*shrink_factor,K)
            # print(floor(heheda)) there is a bug for R package itselt. To be more specific, it is the floor function.
            p_reduce = floor(c(p[1:(d - 1)] * shrink_factor, K))
            # if(p_reduce[1]<=1||p_reduce[2]<=1||p_reduce[3]<=1){ p_reduce=floor(c(p[1:(d-1)]*shrink_factor+1,K)) }
            # if(d==4){ if(p_reduce[3]<=2){ p_reduce[3]=p[3] } if(p_reduce[2]<=2){ p_reduce[2]=3 } if(p_reduce[1]<=2){
            # p_reduce[1]=3 } } if(d==5){ for(test in 1:4) if(p_reduce[test]<=2){ p_reduce[test]=3 } }
            for (test in 1:(d - 1)) {
                if (p_reduce[test] <= 5) {
                  p_reduce[test] = min(p[test], 3)
                }
            }






            # print(p_reduce)
            X_sample_reduce = resize_array(X_sample, c(p_reduce, n))
            r_reduce = r
            # initial parameter downsize
            if (is.na(sum(B0_list[[1]][[1]])) != 1) {
                Number_ini = length(B0_list)
                B0_list_reduce = B0_list
                for (rep in 1:Number_ini) {
                  for (sizepointer in 1:d) {
                    middle = B0_list[[rep]][[sizepointer]]
                    print(dim(middle))
                    print(p_reduce[sizepointer])
                    B0_list_reduce[[rep]][[sizepointer]] = middle[1:p_reduce[sizepointer], 1:r]
                  }
                }
            }


            if (is.na(sum(B0[[1]])) != 1) {
                B0_reduce = B0
                for (sizepointer in 1:d) {
                  B0_reduce = B0[[sizepointer]][1:p_reduce[sizepointer], 1:r]
                }
            }

            if (is.na(sum(BB0[[1]])) != 1) {
                BB0_reduce = resize_array(BB0, c(p_reduce))
            }
            if (is.na(sum(BB0_list[[1]][[1]])) != 1) {
                Number_ini = length(BB0_list)
                BB0_list_reduce = BB0_list
                for (rep in 1:Number_ini) {
                  BB0_list_reduce[[rep]] = resize_array(BB0_list[[rep]], c(p_reduce))
                }

            }

        } else {
            # X_sample_reduce=X_sample
            r_reduce = r
        }
    }


    res = list(X_sample_reduce = X_sample_reduce, r_reduce = r_reduce, B0_reduce = B0_reduce, B0_list_reduce = B0_list_reduce,
        BB0_list_reduce = BB0_list_reduce, BB0_reduce = BB0_reduce, shrink_factor = shrink_factor)
    return(res)
}

uparraysize = function(Xarray_new, pobjectup, index) {
  Xarray = array(0, pobjectup)
  d = length(pobjectup)

  if (d == 2) {
    Xarray[c(index[[1]]), c(index[[2]])] = Xarray_new
  }
  if (d == 3) {
    Xarray[c(index[[1]]), c(index[[2]]), c(index[[3]])] = Xarray_new
  }
  if (d == 4) {
    Xarray[c(index[[1]]), c(index[[2]]), c(index[[3]]), c(index[[4]])] = Xarray_new
  }
  return(Xarray)
}

warmstartupsize = function(warmstart = 3, shrink_factor = 1, p = 1, X_sample_reduce_index = NA, BBreduce = NA) {
    if (warmstart == 1) {
        if (shrink_factor < 1) {
            BBup = uparraysize(Xarray_new = BBreduce, pobjectup = p, index = X_sample_reduce_index[-length(X_sample_reduce_index)])
        } else {
            BBup = BBreduce
        }
    }

    if (warmstart == 3) {
        if (shrink_factor < 1) {
            BBup = resize_array(BBreduce, p)
        } else {
            BBup = BBreduce
        }
    }


    return(BBup)
}

warmstartprint = function(warmstart = 3) {
    cat("Sequential initial stage... \n")
}


predictvalue_nonlinear = function(b0, BB, X, y, family = "gaussian") {
    d = length(dim(X)) - 1
    if (family == "gaussian") {
        n = length(y)
        BB = rTensor::as.tensor(BB)
        if (class(X)[1] == "Tensor") {
        } else {
            X = rTensor::as.tensor(X)
        }
        yhat = c()
        for (i in 1:n) {
            if (d == 4) {
                yhat[i] = rTensor::innerProd(BB, X[, , , , i]) + b0
            }
            if (d == 3) {
                yhat[i] = rTensor::innerProd(BB, X[, , , i]) + b0
            }
        }
        rm(X)
        gc()


    }


    yy = c()
    yyhat = c()
    for (i in 1:n) {
        if (y[i] >= 0) {
            yy[i] = 1
        } else {
            yy[i] = -1
        }
        if (yhat[i] >= 0) {
            yyhat[i] = 1
        } else {
            yyhat[i] = -1
        }


    }
    true = sum(yyhat == yy)
    return(true)

}

predictvalue_linear <- function(b0, BB, X, y, family = "gaussian") {
    if (family == "gaussian") {
        n = length(y)
        if (class(BB)[1] == "Tensor") {
        } else {
            BB = rTensor::as.tensor(BB)
        }


        if (class(X)[1] == "Tensor") {
        } else {
            X = rTensor::as.tensor(X)
        }

        yhat = c()
        for (i in 1:n) {
            yhat[i] = rTensor::innerProd(BB, X[, , , i]) + b0
        }
        yy = c()
        yyhat = c()
        for (i in 1:n) {
            if (y[i] >= 0) {
                yy[i] = 1
            } else {
                yy[i] = -1
            }
            if (yhat[i] >= 0) {
                yyhat[i] = 1
            } else {
                yyhat[i] = -1
            }


        }
        true = sum(yyhat == yy)
        return(true)
    }
}




# -------- norm tensor calculate \int (x-t1)^m1 (x-t2)^m2 dx
interpoly_int <- function(t1, m1, t2, m2, right_end = 1) {
    if (t1 == t2) {
        return(1/(m1 + m2 + 1) * (right_end - t1)^(m1 + m2 + 1))
    }

    min_index <- which.min(c(t1, t2))
    min_t <- c(t1, t2)[min_index]
    max_t <- c(t1, t2)[-min_index]
    min_m <- c(m1, m2)[min_index]
    max_m <- c(m1, m2)[-min_index]
    m1 <- max_m
    t1 <- max_t
    m2 <- min_m
    t2 <- min_t

    s = 0
    for (i in 0:m2) {
        # s = s + choose(m2, i) * 1/(m1 + m2 + 1 - i) * (right_end - t1)^(m1 + m2 + 1 - i) * (t1 - t2)^i
        s = s + choose(m2, i) * 1/(m1 + 1 + i) * ((right_end - t1)^(m1 + 1 + i)) * ((t1 - t2)^(m2 - i))
    }
    return(s)
}

# calculate \int fhat (entry-wise)
fhat_int <- function(coef, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4) {
    num_knots <- length(knots)
    right_end <- knots[num_knots]
    left_end <- knots[1]
    minus_const <- c(rep(left_end, order - 1), knots[1:(num_knots - 2)])
    power_const <- c(1:(order - 1), rep(order - 1, num_knots - 2))
    K <- num_knots - 2 + order - 1
    s = 0
    for (k in 1:K) {
        s = s + coef[k] * 1/(power_const[k] + 1) * (right_end - minus_const[k])^(power_const[k] + 1)
    }
    return(s)
}

# calculate \int fhat^2 (entry-wize)
fhatsqure_int <- function(coef, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4) {
    num_knots <- length(knots)
    K <- num_knots - 2 + order - 1
    right_end <- knots[num_knots]
    left_end <- knots[1]
    minus_const <- c(rep(left_end, order - 1), knots[1:(num_knots - 2)])
    power_const <- c(1:(order - 1), rep(order - 1, num_knots - 2))
    s = 0
    test = matrix(0, K, K)
    for (i in 1:K) {
        for (j in 1:K) {
            test[i, j] = coef[i] * coef[j] * interpoly_int(minus_const[i], power_const[i], minus_const[j], power_const[j],
                right_end = right_end)
            s = s + coef[i] * coef[j] * interpoly_int(minus_const[i], power_const[i], minus_const[j], power_const[j],
                right_end = right_end)
        }
    }
    return(s)
}

# calculate the norm tensor
fhatnorm_ten <- function(BB, knots = c(0, 0.25, 0.5, 0.75, 1), order = 4) {
    pk = dim(BB)
    fhatsqure_intten <- apply(BB, c(1:(length(pk) - 1)), fhatsqure_int, knots = knots, order = order)
    fhat_intten <- apply(BB, c(1:(length(pk) - 1)), fhat_int, knots = knots, order = order)
    norm_ten <- sqrt(fhatsqure_intten - 2 * fhat_intten^2 + (knots[length(knots)] - knots[1]) * fhat_intten^2)
    return(norm_ten)
}

# generate one truncated Spline (not use this due to speed)
truncatS <- function(x, knots, order) {
    part1 <- unlist(lapply(1:(order - 1), function(i) x^i))
    part2 <- unlist(lapply(1:(length(knots) - 2), function(i) ((x - knots[i]) * ((x - knots[i]) > 0))^(order -
        1)))
    return(c(part1, part2))
}






# contract inner product
ctprod <- function(A,B,K){
  da <- dim(A)
  la <- length(da)
  db <- dim(B)
  lb <- length(db)
  Amat <- array(A, dim=c(prod(da[1:(la-K)]),prod(da[(la-K+1):la])))
  Bmat <- array(B, dim=c(prod(db[1:K]),prod(db[(K+1):lb])))
  Cmat <- Amat %*% Bmat
  C <- array(Cmat,dim=c(da[1:(la-K)],db[(K+1):lb]))
  return(C)
}


