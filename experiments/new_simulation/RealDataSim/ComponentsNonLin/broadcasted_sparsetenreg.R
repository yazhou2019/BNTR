
# local dependency
#source("functions_needed.R")
#source("sequential_warmstart.R")
#' @name broadcasted_sparsetenreg
#' @title Fit a broadcasted nonparametric tensor regression model with elastic net regularization
#' @description  A scale-adjusted block-wise descent algorithm to solve the optimization based on Broadcasted Nonparametric Tensor Regression model
#'
#' @param X_sample The training input with basis transformation.
#' @param y The training output.
#' @param r CP rank used in the algorithm
#' @param lambda tuning parameters: \eqn{\lambda_1} in the paper
#' @param alpha tuning parameters: \eqn{\lambda_2} in the paper
#' @param BB0  initial coefficient (tensor form) for the covariate tensor
#' @param beta0 initial intercept term
#' @param Replicates number of initial points used by the algorithm
#' @param restriction used in the future
#' @param knots used in the future
#' @param penalty used in the future
#' @param gamma the tuning parameter for the basis components in CP decomposition (equivalent to \eqn{\lambda_1} of the elastic net)
#' @param alpha_gamma the tuning parameter for the basis components in CP decomposition (equivalent to \eqn{\lambda_2} of the elastic net)
#' @param family used in the future
#' @param Z_sample used in the future
#' @param B0 ipitial coefficient (CP decomposition form) for the covariate tensor
#' @param epsilon tolerance value in block update
#' @param MaxIter maximum of the outer iteration number
#' @param TolFun tolerance value of the outer iteration (the whole algorithm)
#' @param Memorysave run the algorithm with saving (1) or withougt saving the memory (0)
#' @param rescale_l2 initial rescale strategy: l_2 method (1) or l_1 method (0).
#' @param rescale_all rescale all CP components including the coefficient of the basis (1) or without the coeffcients of basis (0)
#' @param manifold used in the future
#' @param QCQP update the basis components by QCQP optimization (1) or not (0)
#' @param warmstart warmstart strategy; 1 is sequential; 3 is one time.
#' @param improvement numerical stablization improvementation (1) or not (0)
#' @param stoppingrule_penalty stoppoing rule. 1 means that the objective function is loss + penalty; 0 means that the objective function is just loss.
#' @param  shrink_factor_number number of sequential initialization
#' @param startsize the minimal down-size in the sequential initialization
#' @return
#' \item{beta0}{the estimated intercept term}
#' \item{beta}{the estimated coefficient tensor (in CP decomposition form)}
#' \item{ob_f_loss_penalty_final}{the final objective function value track}
#' \item{others}{used for the package developer}
#' @seealso validation_broadcasted_sparsetenreg, sequential_warmstart
#' @references Y. Zhou, R. K. W. Wong and K. He. Broadcasted Nonparametric Tensor Regression
#' @examples
#' data("X_horse")
#' BB <- X_horse
#'
#' set.seed(2019)
#'
#' # sample size of the training and test set
#' n_train <- 400
#' n_test <- 100
#'
#' # signal level
#' signal_level = 0.5
#'
#'
#' # the tuning parameters have been tuned by authors though validation; see the paper.
#' # rank used in the algorithm
#' rank <- 4
#' # the tuning parameters of the elastic net
#' lambda_1 <- 5
#' lambda_2 <- 0.5  # 0 corresponds to lasso; 1 corresponds to ridge
#' # use the initial point tuned based on 400 training, or not
#' input_initial <- 0
#'
#'
#' # training data
#' X_train = array(runif(prod(c(dim(BB), n_train)), 0, 1), c(dim(BB), n_train))
#' # broadcated procedure
#' BroX_train = X_train + 0.6 * sin(2 * pi * (X_train - 0.5)^2)
#'
#' y_train = 1 + crossprod(matrix(BroX_train, c(prod(dim(BB)), n_train)), as.vector(BB)) + signal_level * rnorm(n_train)
#'
#' # Transform to truncated power basis
#' num_knots <- 5
#' knots = quantile(c(X_train), probs = c(seq(0, 1, 1/(num_knots - 1))))
#' tildePhiX_train = tildePhiX_trans(X_train, knots)
#'
#' # BroadcasTR
#' if (input_initial == 1) {
#'  data("initial_point")
#'  beta0 <- initial_point$beta0
#'  B0 <- initial_point$B0
#'  warmstart <- 0
#' } else {
#' beta0 <- NA
#'  B0 <- NA
#'  warmstart <- 1
#' }
#' res = broadcasted_sparsetenreg(tildePhiX_train, y_train, r = rank, lambda = lambda_1, alpha = lambda_2, warmstart = warmstart, beta0 = beta0, B0=B0, Replicates=1)
#'

#' # norm tensor
#' normtensor <- fhatnorm_ten(full_R(res$beta), knots)
#'
#' # plot the true and estimated region
#' par(mar = c(0, 0, 0, 0), mfcol = c(1, 2), mai = c(0.2, 0.4, 0.4, 0.2))
#' plot(c(1, dim(BB)[1]), c(1, dim(BB)[2]), xaxt = "n", yaxt = "n", type = "n")
#' rasterImage(BB, 0, 0, dim(BB)[1] + 1, dim(BB)[2] + 1)
#' mtext("True coefficient tensor (matrix)", 3, line = 0.2)
#' plot(c(1, dim(BB)[1]), c(1, dim(BB)[2]), xaxt = "n", yaxt = "n", type = "n")
#' mtext("Estimated norm tensor (matrix)", 3, line = 0.2)
#' rasterImage(normtensor, 0, 0, dim(BB)[1] + 1, dim(BB)[2] + 1)
#'
#'
#'
#'
#' # generate the test data
#' X_test = array(runif(prod(c(dim(BB), n_test)), 0, 1), c(dim(BB), n_test))
#' BroX_test = X_test + 0.6 * sin(2 * pi * (X_test - 0.5)^2)
#' y_test = 1 + crossprod(matrix(BroX_test, c(prod(dim(BB)), n_test)), as.vector(BB)) + signal_level * rnorm(n_test)
#'
#'
#' # prediction on the test set
#' tildePhiX_test = tildePhiX_trans(X_test, knots)
#' y_pre = res$beta0 + crossprod(matrix(tildePhiX_test, c(prod(dim(full_R(res$beta))), n_test)), as.vector(full_R(res$beta)))
#'
#' # prediction performance
#' cat("The prediction performance in these tunning parameters: \n", "MSPE =", sum((y_test - y_pre)^2)/n_test, "\n")
#'
#' #set.seed(1224) # used for validation method
#' #n_vali = 0.25 * n_train
#' # generate the test data
#' # X_vali = array(runif(prod(c(dim(BB), n_vali)), 0, 1), c(dim(BB), n_test))
#' # BroX_vali = X_vali + 0.6 * sin(2 * pi * (X_vali - 0.5)^2)
#' # y_vali = 1 + crossprod(matrix(BroX_vali, c(prod(dim(BB)), n_vali)), as.vector(BB)) + signal_level * rnorm(n_vali)
#'
#'
#' # R=c(1,2,3,4,5)
#' # alpha=c(0,0.5,1)
#' # lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)
#'
#' # res_valiandtest <- validation_broadcasted_sparsetenreg(R,alpha,lambda,X_train,y_train,X_vali,y_vali,X_test,y_test, num_knots=5, order=4)




broadcasted_sparsetenreg <- function(X_sample, y, r = 2,lambda = 0, alpha = 0.5,
    gamma = 0, alpha_gamma = 0, warmstart = 1, beta0 = NA, B0 = NA, BB0 = NA, Replicates = 1, QCQP = 1,
    epsilon = 1e-07, MaxIter = 10000, TolFun = 1e-04,  rescale_l2 = 1, rescale_all = 0, Memorysave = 0, manifold = 0,
     family = "gaussian", Z_sample = NA, improvement = 1, stoppingrule_penalty = 1, shrink_factor_number = 5, startsize = 3, restriction = 0, knots = 0, penalty = c("L1L2")) {
    # This is the main function of Broadcasted Nonparametric Tensor Regression (BroadcaTR).

    # parse inputs:
    ###############commonly used##########################
    # X_sample: the basis tensor of the original input tensor
    # y: the output
    # r: CP rank used in the algorithm
    # lambda: lambda_1
    # aloha: lambda_2
    # BB0: initial point in the tensor form
    # beta0: initial intercept
    # Replicates: number of initial points used by the algorithm


    ################ parameters can be used by default##########
    # restriction: used in the future
    # knots: used in the future
    # penalty: used in the future
    # gamma: the tuning parameter for the basis components in CP
    # alpha_gamma: the tuning parameter for the basis components in CP
    # family: used in the future
    # Z_sample: used in the future
    # B0:initial point in the CP components form
    # epsilon: tolerance value in glmnet
    # MaxIter: maximum of the outer iteration number
    # TolFun: tolerance value of the whole algorithm
    # Memorysave: run the algorithm with saving or withougt saving the memory
    # rescale_l2: initial rescale strategy
    # rescale_all: rescale all components or just the coeffcients components
    # manifold: update the basis components by Oblique manifold optimization
    # QCQP: update the basis components by QCQP manifold optimization
    # warmstart: warmstart strategy. 1 is sequential; 3 is one time
    # improvement: numerical stablization improvementation
    # stoppingrule_penalty: stoppoing rule. 1 is loss + penalty; 0 is loss.
    # shrink_factor_number: number of sequential initialization
    # startsize: the minimal size in the sequential initialization


    if (alpha != 0 && alpha != 1) {
        elastic_indicator <- 1
    } else {
        elastic_indicator <- 0
    }

    # used to improve some converge problems of glmnet function
    control_iteration <- 0

    # fix the scale of the norm of coefficient of basis to be 1 to match the paper.  it may have nummerical effects
    # in applications
    constantscale <- 1

    # used for the warmstart 3
    shrink_factor_scale <- 5

    # used for feasible adjustment for QCQP
    norm1_iterationstart <- 1


    n = length(y)
    lambda = lambda/n
    p1 = dim(X_sample)
    p = p1[-length(p1)]

    # used to check convergence
    Convergence_Index = 0


    if (is.na(Z_sample) == 1) {
        X = rep(1, n)
    } else {
        zdim = dim(Z_sample)[2]
        X = matrix(0, n, zdim + 1)
        X[, 1:zdim] = Z_sample
        X[, zdim + 1] = rep(n, 1)
    }
    X = as.matrix(X)
    xdim = dim(X)[2]  # the dimension of other parameters

    # transform the input to tensor object
    d = length(p)
    if (class(X_sample)[1] == "Tensor") {
        TM = X_sample
    } else {
        TM = rTensor::as.tensor(X_sample)
    }



    beta = list()
    # one time down-size warmstart
    if (warmstart == 3) {
        fit_singlewarm = single_warmstart(TolFun = 0.01, rescale_all = 1, rescale_l2 = 1, r = r, X_sample = X_sample,
            y = y, lambda = 0, gamma = 0, alpha_gamma = 0, alpha = 0, Replicates = 5, family = family, manifold = 0,
            warmstart = 3, QCQP = 0, shrink_factor_scale = shrink_factor_scale, Memorysave = 0)
        beta = fit_singlewarm$upsizebeta
        beta0 = fit_singlewarm$bbeta0
    }

    # sequential down-size warmtsart.
    if (warmstart == 1) {
        cat(rep("+", 30), " \n")
        cat(rep("+", 30), " \n")
        cat(rep("+", 6), "START::sequential initial procedure", rep("+", 6), "\n")
        fit_seqwarm = sequential_warmstart(TolFun = 0.01, r = r, X_sample = X_sample, y = y, penalty = c("L1"),
            lambda = 0, gamma = 0, alpha_gamma = 0, alpha = 0, family = "gaussian", shrink_factor_number = shrink_factor_number,
            startsize = startsize, Memorysave = Memorysave, MaxIter=50)
        beta = fit_seqwarm$bbeta
        beta0 = fit_seqwarm$bbeta0
        cat(rep("-", 6), "END::sequential initial procedure", rep("-", 7), "\n")
        cat(rep("-", 30), " \n")
        cat(rep("-", 30), " \n")
        cat(rep("+", 30), " \n")
        cat(rep("+", 30), " \n")
        cat(rep("+", 10), "START::the algorithm", rep("+", 10), "\n")
    }

    rm(X_sample)
    gc()

    # speed up if memory is large enough; speed up by sacrificing the memory
    Md = list()
    if (Memorysave != 1) {
        for (dd in 1:d) {
            Md[[dd]] = rTensor::unfold(TM, c(d + 1, dd), c(1:d)[-dd])
        }
        rm(TM)
        gc()
    }

    # final loss initialization
    dev_final = Inf
    # finial loss + penalty initialization
    dev_final_penalty = Inf
    # Objective function: just loss
    ob_f_final = 0
    # # Objective function: loss + penalty
    ob_f_loss_penalty_final = 0

    # loss of different initial points
    ob_f_differentinitial = c()
    # loss + penalty of different initial points (may be used in the future)
    ob_f_loss_penalty_differentinitial = c()
    # linear predictor of different initial points
    linear_predictor_differentinitial = list()

    # The iteration number
    num_final = 0

    # If Replicates > 1, then the will be Replicates - 1 random initial points
    for (rep in 1:Replicates) {

        # Objective function final

        ob_f = c()
        ob_f_loss_penalty = c()
        # initialize the CP components from uniform [-1,1] or standard norm distribution
        set.seed(rep * 10)
        # when use warmstart, we also can use random initial points and compare the final objective function value, and
        # choose the best.
        if (rep > 1) {
            B0 = NA
            BB0 = NA
            warmstart = 0
        }



        if (is.na(sum(B0[[1]])) == 1 & is.na(sum(BB0[[1]])) == 1) {
            if (warmstart != 0) {
                if (QCQP == 1 || manifold == 1) {
                  # to fix the norm 1
                  for (i in 1:r) {
                    norm1scale = sqrt(sum(beta[[d]][, i]^2)) * constantscale
                    beta[[d]][, i] = beta[[d]][, i]/norm1scale
                    beta[[1]][, i] = beta[[1]][, i] * norm1scale
                  }
                }

            } else {
                for (i in 1:d) {
                  # both partern is OK
                  beta[[i]] = matrix(stats::runif(p[i] * r, -1, 1), p[i], r)
                  # beta[[i]]=matrix(rnorm(p[i]*r),p[i],r)
                }
                if (QCQP == 1 || manifold == 1) {
                  # to fix the norm 1
                  for (i in 1:r) {
                    norm1scale = sqrt(sum(beta[[d]][, i]^2)) * constantscale
                    beta[[d]][, i] = beta[[d]][, i]/norm1scale
                    beta[[1]][, i] = beta[[1]][, i] * norm1scale
                  }
                }
            }
        } else {
            if (warmstart != 0) {
                # If there is a warmstart, then just use the warmstart initial point
            } else {
                if (is.na(sum(B0[[1]])) != 1) {
                  for (i in 1:d) {
                    beta[[i]] = matrix(B0[[i]], p[i], r)
                  }
                }
                if (is.na(sum(BB0[[1]])) != 1) {
                  cpBB0 = rTensor::cp(rTensor::as.tensor(BB0), r)
                  beta = cpBB0$U
                  cp_lambda = cpBB0$lambdas
                  for (i in 1:d) {
                    for (rr in 1:r) {
                      betaii = beta[[i]]
                      beta[[i]][, rr] = betaii[, rr] * (cp_lambda[rr])^(1/d)
                    }
                  }

                }



            }

            if (QCQP == 1 || manifold == 1) {
                # fix the norm to be 1 to match the paper
                for (i in 1:r) {
                  norm1scale = sqrt(sum(beta[[d]][, i]^2)) * constantscale
                  beta[[d]][, i] = beta[[d]][, i]/norm1scale
                  beta[[1]][, i] = beta[[1]][, i] * norm1scale
                }
            }

        }




        # main loop
        for (iter in 1:MaxIter) {
            num = iter
            if (iter == 1) {
                # initial value of the intercept (if without warmstart or user's input)
                if (is.na(sum(beta0)) == 1) {
                  if (warmstart == 0) {
                    df <- data.frame(X, y)
                    fit <- stats::glm(y ~ . - 1, family = family, data = df, control = list(maxit = 1e+06, epsilon = 1e-08))
                    beta0 = stats::coef(fit)
                  }
                } else {
                  beta0 = beta0
                }



                linear_predictor = as.matrix(X) %*% beta0
                dev0 = loglike(y, linear_predictor, family = family)
                dev0_penalty = dev0 + n * lambda * penaltyvalue(beta = beta, alpha = alpha, manifold = manifold)
                ob_f[1] = 0
                ob_f_loss_penalty[1] = 0

            } else {
                # stopping rule is in this part
                eta = Xj %*% as.vector(beta[[d]])
                X_eta = matrix(0, n, xdim + 1)
                X_eta[, 1:xdim] = X
                X_eta[, xdim + 1] = eta
                df <- data.frame(X_eta, y)

                # used to avoid the numerical problems
                if (sum(abs(eta)) < 1e-16) {
                  message("Warning: beyond the numerical boundary")
                  training_error = loglike(y, linear_predictor, family = family)
                  res = list(beta0 = beta0, beta = beta, dev = dev0, error = training_error, num = num)
                  return(res)
                }

                if (improvement == 1 && iter > control_iteration) {
                  df <- data.frame(X, y)
                  fit <- stats::glm(y ~ . - 1, family = family, offset = eta, data = df, control = list(maxit = 1e+07,
                    epsilon = epsilon))
                  beta_eta = c(stats::coef(fit), 1)
                  end = length(beta_eta)
                  beta0 = beta_eta[-end]
                } else {
                  fit <- stats::glm(y ~ . - 1, family = family, data = df, control = list(maxit = 1e+07, epsilon = epsilon))
                  beta_eta = stats::coef(fit)
                  end = length(beta_eta)
                  beta0 = beta_eta[-end]
                }








                # numerical stabilization improvement
                if (improvement == 1 && iter > control_iteration) {
                  linear_predictor = eta + beta0
                  devtmp = loglike(y, linear_predictor, family = family)

                } else {
                  linear_predictor = X_eta %*% beta_eta
                  devtmp = loglike(y, linear_predictor, family = family)
                }


                if (improvement == 1 && iter > control_iteration) {
                  # when updating the intercept, do not update the scale of CP components.
                  if (iter > norm1_iterationstart) {
                    # feasible adjustment and fix norm 1 if it is infeasible, then it will be feasible if it if feasible,the the
                    # objective funtion (with penalty) will be smaller
                    if (QCQP == 1 || manifold == 1) {
                      # to fix the norm 1
                      for (i in 1:r) {
                        norm1scale = sqrt(sum(beta[[d]][, i]^2)) * sqrt(constantscale)
                        if (norm1scale > 0.01) {
                          beta[[d]][, i] = beta[[d]][, i]/norm1scale
                          beta[[1]][, i] = beta[[1]][, i] * norm1scale
                        }
                      }
                    }
                  }

                  beta = rescale(beta, rescale_all = rescale_all, rescale_l2 = rescale_l2, elastic = elastic_indicator,
                    alpha = alpha)
                } else {
                  # when updating the intercept,update the scale of CP components for numerical stabalization.
                  beta = arrange(beta, beta_eta[length(beta_eta)])
                  beta = rescale(beta, rescale_all = rescale_all, rescale_l2 = rescale_l2, elastic = elastic_indicator,
                    alpha = alpha)
                }


                ob_f[(d + 1) * (iter - 1) + 1] = devtmp
                ob_f_loss_penalty[(d + 1) * (iter - 1) + 1] = devtmp + n * lambda * penaltyvalue(beta = beta, alpha = alpha,
                  manifold = manifold, QCQP = QCQP)
                devtmp_penalty = devtmp + n * lambda * penaltyvalue(beta = beta, alpha = alpha, manifold = manifold,
                  QCQP = QCQP)


                ###### used to numerical boundary checking. whether the linear_predictor is too big
                if (is.nan(devtmp) == 1 || is.na(devtmp) == 1 || is.na(devtmp_penalty) == 1 || is.nan(devtmp_penalty) ==
                  1) {
                  # In the former iteration, it has been reach the boundary
                  message("Warning: beyond the numerical boundary")
                  training_error = loglike(y, linear_predictor, family = family)
                  res = list(beta0 = beta0, beta = beta, dev = dev0, error = training_error, num = num, ob_f_final = ob_f/n,
                    ob_f_loss_penalty_final = ob_f_loss_penalty/n)
                  return(res)
                }

                diffdev = devtmp - dev0
                dev0 = devtmp

                # Two stopping rule, loss or loss + penalty
                diffdev_penalty = devtmp_penalty - dev0_penalty
                dev0_penalty = devtmp_penalty
                if (stoppingrule_penalty == 1) {
                  v1 = abs(diffdev_penalty)
                  v2 = TolFun * (abs(dev0_penalty) + 1)
                } else {
                  v1 = abs(diffdev)
                  v2 = TolFun * (abs(dev0) + 1)
                }



                if (v1 < v2) {
                  Convergence_Index = 1

                  if(improvement==1&&iter> control_iteration){
                    beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)
                  }else{
                    beta=arrange(beta,beta_eta[length(beta_eta)])
                    beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)
                  }

                  cat("Congratulation!!!reach the stopping rule \n")
                  break
                }

                # when the iteration number is big, just update the intercept; while not big, for numerical reason, update them
                # togeter
                if (improvement == 1 && iter > control_iteration) {
                  beta = rescale(beta, rescale_all = rescale_all, rescale_l2 = rescale_l2, elastic = elastic_indicator,
                                 alpha = alpha)
                } else {
                  beta = arrange(beta, beta_eta[length(beta_eta)])
                  beta = rescale(beta, rescale_all = rescale_all, rescale_l2 = rescale_l2, elastic = elastic_indicator,
                                 alpha = alpha)
                }


            }

            # Alternatively update the block
            eta0 = as.matrix(X) %*% beta0
            # loop for CP components
            for (j in 1:d) {
                if (j == 1) {
                  cumkr = t(as.matrix(rep(1, r)))
                }
                if (j == d) {
                  # save memory method
                  if (Memorysave == 1) {
                    Md[[j]] = rTensor::unfold(TM, c(d + 1, j), c(1:d)[-j])
                    Xj = matrix(Md[[j]]@data %*% cumkr, n, p[j] * r)
                    Md[[j]] = 0
                    gc()
                  } else {
                    Xj = matrix(Md[[j]]@data %*% cumkr, n, p[j] * r)
                  }
                } else {
                  if (j == (d - 1)) {
                    if (Memorysave == 1) {
                      Md[[j]] = rTensor::unfold(TM, c(d + 1, j), c(1:d)[-j])
                      Xj = matrix(Md[[j]]@data %*% rTensor::khatri_rao(beta[[d]], cumkr), n, p[j] * r)
                      Md[[j]] = 0
                      gc()
                    } else {
                      Xj = matrix(Md[[j]]@data %*% rTensor::khatri_rao(beta[[d]], cumkr), n, p[j] * r)
                    }
                  } else {
                    if (Memorysave == 1) {
                      Md[[j]] = rTensor::unfold(TM, c(d + 1, j), c(1:d)[-j])
                      Xj = matrix(Md[[j]]@data %*% rTensor::khatri_rao(rTensor::khatri_rao_list(beta[c(d:(j + 1))]),
                        cumkr), n, p[j] * r)
                      Md[[j]] = 0
                      gc()
                    } else {
                      Xj = matrix(Md[[j]]@data %*% rTensor::khatri_rao(rTensor::khatri_rao_list(beta[c(d:(j + 1))]),
                        cumkr), n, p[j] * r)
                    }
                  }
                }

                Xj_eta = matrix(0, n, p[j] * r + 1)
                # Corresponding to the CP components
                Xj_eta[, 1:(p[j] * r)] = Xj
                # Corresponding to the intercept
                Xj_eta[, p[j] * r + 1] = eta0
                # print(Xj_eta) used for boundary beyound, Xj
                if (sum(is.nan(Xj)) != 0) {
                  # In the former iteration, it has been reach the boundary
                  message("Warning: beyond the number boundary in the software")
                  training_error = loglike(y, linear_predictor, family = family)
                  res = list(beta0 = beta0, beta = beta, dev = dev0, error = training_error, num = num, ob_f_final = ob_f/n,
                    ob_f_loss_penalty_final = ob_f_loss_penalty/n)
                  return(res)
                }

                if (sum(Xj) == 0) {
                  # In the former iteration, it has been reach the boundary
                  message("Warning: lambda is too big")
                  training_error = loglike(y, linear_predictor, family = family)
                  res = list(beta0 = beta0, beta = beta, dev = dev0, error = training_error, num = num, ob_f_final = ob_f/n,
                    ob_f_loss_penalty_final = ob_f_loss_penalty/n)
                  return(res)
                }








                if (j < d) {
                  # below: for L1 or L2 penalty
                  {


                    if (improvement == 1 && iter > control_iteration) {

                      fit <- glmnet::glmnet(Xj_eta[, -(p[j] * r + 1)], y - Xj_eta[, (p[j] * r + 1)], family = family,
                        lambda = lambda, intercept = 1, alpha = alpha, maxit = 1e+08, standardize = FALSE, thresh = epsilon)
                      Bd_new_vec_and = glmnet::coef.glmnet(fit, exact = TRUE)[-1]



                      beta[[j]] = matrix(Bd_new_vec_and, p[j], r)
                      eta0 = Xj_eta[, (p[j] * r + 1)] + glmnet::coef.glmnet(fit, exact = TRUE)[1]

                      linear_predictor = Xj_eta[, -(p[j] * r + 1)] %*% Bd_new_vec_and + eta0
                      devtmp = loglike(y, linear_predictor, family = family)


                    } else {
                      # This is the current function a test
                      p.fac = rep(1, p[j] + 1)
                      p.fac[p[j] + 1] = 0

                      # intercept is in (p[j]+1)-th colume
                      fit <- glmnet::glmnet(Xj_eta, y, family = family, lambda = lambda, intercept = 0, penalty.factor = p.fac,
                        alpha = alpha, maxit = 1e+08, standardize = FALSE, thresh = epsilon)
                      Bd_new_vec_and = glmnet::coef.glmnet(fit, exact = TRUE)[-1]

                      # lambda_1 is too big to push the components to be zero
                      if (sum(abs(Bd_new_vec_and[-length(Bd_new_vec_and)])) == 0) {
                        message("Warning: lambda is too big, Bj, j<d, return 0")
                        training_error = loglike(y, linear_predictor, family = family)
                        res = list(beta0 = beta0, beta = beta, dev = dev0, error = training_error, ob_f_final = ob_f/n,
                          ob_f_loss_penalty_final = ob_f_loss_penalty/n)
                        return(res)
                      }

                      beta[[j]] = matrix(Bd_new_vec_and[-length(Bd_new_vec_and)], p[j], r)
                      eta0 = eta0 * Bd_new_vec_and[length(Bd_new_vec_and)]

                      # The objective values
                      linear_predictor = Xj_eta %*% Bd_new_vec_and + eta0
                      devtmp = loglike(y, linear_predictor, family = family)


                    }

                  }
                  # above: for L1 or L2 penalty

                  ob_f[(d + 1) * (iter - 1) + 1 + j] = devtmp
                  devtmp_penalty = devtmp + n * lambda * penaltyvalue(beta = beta, alpha = alpha, manifold = manifold,
                    QCQP = QCQP)
                  ob_f_loss_penalty[(d + 1) * (iter - 1) + 1 + j] = devtmp_penalty

                } else {

                  if (QCQP == 1) {
                    # QCQP for the d-th mode
                    betajstandard = beta[[j]]
                    linearpartjstandard = Xj_eta[, -(p[j] * r + 1)] %*% c(beta[[j]]) + Xj_eta[1, (p[j] * r + 1)]
                    devtmpstandard = loglike(y, linearpartjstandard, family = family)
                    fit = dualascentforalpha(y = y - Xj_eta[, (p[j] * r + 1)], DX = Xj_eta[, -(p[j] * r + 1)],
                      alphaold = beta[[j]], C = constantscale, Maxiterfordual = 1000, TolFundualascent = 0.001,
                      MaxIterInner = 40)
                    beta[[j]] = matrix(fit$alphavec_final, p[j], r)
                    linear_predictor = Xj_eta[, -(p[j] * r + 1)] %*% c(beta[[j]]) + Xj_eta[1, (p[j] * r + 1)]
                    devtmp = loglike(y, linear_predictor, family = family)

                    if (devtmp > devtmpstandard) {
                      beta[[j]] = betajstandard
                      linear_predictor = linearpartjstandard
                      devtmp = devtmpstandard
                    }




                  } else {
                    # penalize the d-th mode
                    if (improvement == 1 && iter > control_iteration) {
                      fit <- glmnet::glmnet(Xj_eta[, -(p[j] * r + 1)], y - Xj_eta[, (p[j] * r + 1)], family = family,
                        lambda = gamma, intercept = 1, alpha = alpha_gamma, maxit = 1e+07, standardize = FALSE,
                        thresh = epsilon)
                      Bd_new_vec_and = glmnet::coef.glmnet(fit, exact = TRUE)[-1]

                      beta[[j]] = matrix(Bd_new_vec_and, p[j], r)
                      eta0 = Xj_eta[, (p[j] * r + 1)]

                      linearpart = Xj_eta[, -(p[j] * r + 1)] %*% Bd_new_vec_and + eta0 + glmnet::coef.glmnet(fit, exact = TRUE)[1]
                      devtmp = loglike(y, linearpart, family = family)

                    } else {
                      p.fac = rep(1, p[j] + 1)
                      p.fac[p[j] + 1] = 0

                      fit <- glmnet::glmnet(Xj_eta, y, family = family, lambda = gamma, intercept = 0, penalty.factor = p.fac,
                        alpha = alpha_gamma, maxit = 1e+07, standardize = FALSE, thresh = epsilon)
                      Bd_new_vec_and = glmnet::coef.glmnet(fit, exact = TRUE)[-1]


                      beta[[j]] = matrix(Bd_new_vec_and[-length(Bd_new_vec_and)], p[j], r)
                      eta0 = eta0 * Bd_new_vec_and[length(Bd_new_vec_and)]
                      linearpart = Xj_eta %*% Bd_new_vec_and + eta0
                      devtmp = loglike(y, linearpart, family = family)
                    }

                  }

                  # the objective values
                  ob_f[(d + 1) * (iter - 1) + 1 + j] = devtmp
                  devtmp_penalty = devtmp + n * lambda * penaltyvalue(beta, alpha, manifold = manifold, QCQP = QCQP)
                  ob_f_loss_penalty[(d + 1) * (iter - 1) + 1 + j] = devtmp_penalty



                }


                cumkr = rTensor::khatri_rao(beta[[j]], cumkr)


            }


        }

        # the final linear predictor of different initial points
        linear_predictor_differentinitial[[rep]] = linear_predictor

        # the final loss of different initial points.
        ob_f_differentinitial[rep] = dev0
        ob_f_loss_penalty_differentinitial[rep] = dev0_penalty

        # used for the stopping rule
        if (stoppingrule_penalty == 1) {
            # the loss+penalty
            v3 = dev0_penalty
            v4 = dev_final_penalty
        } else {
            # just the loss, used by TensorReg
            v3 = dev0
            v4 = dev_final
        }



        # record if it has a smaller objective value
        if (v3 < v4) {
            beta0_final = beta0
            beta_final = beta
            dev_final = dev0
            linear_predictor_final = linear_predictor
            dev_final_penalty = ob_f_loss_penalty[length(ob_f_loss_penalty)]
            # final iteration number for a specific initial point.
            num_final = num

            # get the final objective funtion that we used rescale by the sample size n
            ob_f_final = ob_f/n
            ob_f_loss_penalty_final = ob_f_loss_penalty/n

        }

    }


    training_error = loglike(y, linear_predictor_final, family = family)



    if (Convergence_Index != 1)
        message("do NOT reach the stopping rule; need more iteration, or smaller tolerance")




    res = list(beta0 = beta0_final, beta = beta_final, dev = dev_final, error = training_error, num = num_final,
        ob_f_final = ob_f_final, ob_f_loss_penalty_final = ob_f_loss_penalty_final, ob_f_differentinitial = ob_f_differentinitial,
        ob_f_loss_penalty_differentinitial = ob_f_loss_penalty_differentinitial, linear_predictor_differentinitial = linear_predictor_differentinitial)
    return(res)
}













