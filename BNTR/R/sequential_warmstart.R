#' @title sequential_warmstart
#' @description  sequential warmstart strategies
#' @seealso validation_broadcasted_sparsetenreg, broadcasted_sparsetenreg
#' @references Section D of the supplementary for Y. Zhou, R. K. W. Wong and K. He. Broadcasted Nonparametric Tensor Regression

sequential_warmstart = function(X_sample, y, r = 2, restriction = 0, knots = 0, penalty = c("L1L2"), lambda = 0,
    gamma = 0, alpha_gamma = 0, alpha = 0, family = "gaussian", Z_sample = NA, B0 = NA, BB0 = NA, Replicates = 5,
    epsilon = 1e-07, MaxIter = 20, TolFun = 0.01, Memorysave = 0, rescale_l2 = 1, rescale_all = 0, manifold = 0,
    improvement = 1, stoppingrule_penalty = 1, QCQP = 1, beta0 = NA, beta0_vec = NA, shrink_factor_number = 2,
    startsize = 3) {
    # startsize is the minimal downsize shrink_factor_number is the number of different downsize

    fitlist_warmstart = list()
    p = dim(X_sample)
    n = p[length(p)]
    p = p[-length(p)]
    common_diff = (p[1] - startsize)/(shrink_factor_number - 1)
    shrink_factor_seq = seq(startsize, p[1], common_diff)/p[1]
    shrink_factor_scale_seq = n/(r * sum(p) * shrink_factor_seq)
    num_sequential = length(shrink_factor_scale_seq)

    for (iter in 1:num_sequential) {
        if (iter == 1) {
            BB0_initial = NA
            beta0_initial = NA
        } else {
            BB0_initial = full_R(bbeta)
            beta0_initial = bbeta0
        }
        cat(stringr::str_c(iter, "-th initialization in the sequential procedure"), "\n")
        fitlist_warmstart[[iter]] = single_warmstart(TolFun = TolFun, rescale_all = rescale_all, rescale_l2 = rescale_l2,
            r = r, X_sample = X_sample, y = y, lambda = lambda, gamma = gamma, alpha_gamma = alpha_gamma, alpha = alpha,
            Replicates = 1, family = family, manifold = 0, warmstart = 3, QCQP = 0, BB0 = BB0_initial, beta0 = beta0_initial,
            shrink_factor_scale = shrink_factor_scale_seq[iter], Memorysave = Memorysave)
        bbeta = fitlist_warmstart[[iter]]$upsizebeta
        bbeta0 = fitlist_warmstart[[iter]]$downsizeres$beta0



    }
    res = list(bbeta = bbeta, bbeta0 = bbeta0)
    return(res)
}

#' @title single_warmstart
#' @description  one times down-size warmstart strategies, which will be used in the sequential warmstart.
#' @seealso validation_broadcasted_sparsetenreg, broadcasted_sparsetenreg
#' @references Section D of the supplementary for Y. Zhou, R. K. W. Wong and K. He. Broadcasted Nonparametric Tensor Regression

single_warmstart <- function(X_sample, y, r = 2, restriction = 0, knots = 0, penalty = c("L1L2"), lambda = 0, gamma = 0,
    alpha_gamma = 0, alpha = 0, family = "gaussian", Z_sample = NA, B0 = NA, BB0 = NA, Replicates = 5, epsilon = 1e-07,
    MaxIter = 1000, TolFun = 0.01, Memorysave = 0, rescale_l2 = 1, rescale_all = 0, manifold = 0, warmstart = 3,
    improvement = 1, stoppingrule_penalty = 1, QCQP = 1, beta0 = NA, beta0_vec = NA, shrink_factor_scale = 5) {


    # MaxIter=1000 TolFun=0.0001
    n = length(y)
    p1 = dim(X_sample)
    p = p1[-length(p1)]




    B0_reduce = NA
    BB0_reduce = NA



    if (is.na(Z_sample) == 1) {
        X = rep(1, n)
    } else {
        z2 = dim(Z_sample)[2]
        X = matrix(0, n, z2 + 1)
        X[, 1:z2] = Z_sample
        X[, z2 + 1] = rep(n, 1)
    }
    X = as.matrix(X)
    x2 = dim(X)[2]

    d = length(p)
    if (class(X_sample)[1] == "Tensor")
        TM = X_sample else TM = rTensor::as.tensor(X_sample)

    p_reduce = p

    # In the future we can consider more initial strategies.
    if (warmstart != 0) {
        K = p[length(p)]
        # this warmstart method is one times down size
        {


            wanmstartdowsize_res = warmstartdowsize(warmstart = warmstart, r = r, X_sample = X_sample, shrink_factor_scale = shrink_factor_scale,
                B0 = B0, BB0 = BB0)

            X_sample_reduce = wanmstartdowsize_res$X_sample_reduce
            r_reduce = wanmstartdowsize_res$r_reduce
            B0_reduce = wanmstartdowsize_res$B0_reduce

            BB0_reduce = wanmstartdowsize_res$BB0_reduce

            shrink_factor = wanmstartdowsize_res$shrink_factor


            ### where the adjustment

            BBreduce = broadcasted_sparsetenreg(X_sample_reduce, y, r = r_reduce, penalty = c("L1L2"), lambda = 0,
                gamma = 0, alpha_gamma = 0, alpha = 0, family = "gaussian", Z_sample = NA, Replicates = Replicates,
                epsilon = 1e-07, MaxIter = 1000, TolFun = TolFun, Memorysave = Memorysave, rescale_l2 = 1, rescale_all = 1,
                manifold = 0, QCQP = 0, B0 = B0_reduce, BB0 = BB0_reduce, beta0 = beta0, warmstart = 0)

            testlist = BBreduce

            # return(BBreduce)
            beta0 = BBreduce$beta0
            BBreduce = full_R(BBreduce$beta)
            BBup = warmstartupsize(warmstart = warmstart, shrink_factor = shrink_factor, p = p,
                BBreduce = BBreduce)



        }

        # warmstartprint(warmstart = 3)

        cpD = rTensor::cp(rTensor::as.tensor(BBup), r_reduce)
        # return(cpD)
        beta = cpD$U
        cp_lambda = cpD$lambdas
        # return(beta)
        for (i in 1:d) {
            for (rr in 1:r_reduce) {
                betaii = beta[[i]]
                beta[[i]][, rr] = betaii[, rr] * (cp_lambda[rr])^(1/d)
            }
        }




    } else {
        beta = list()
    }


    test_list_return = list(downsizeres = testlist, upsizebeta = beta, shrink_factor = shrink_factor, bbeta0 = beta0)
    return(test_list_return)

}




