#' @title Tune parameters through validation method
#' @name validation_broadcasted_sparsetenreg
#' @param R grids of the rank used in the algorithm
#' @param alpha grids of \eqn{\lambda_2} used in the algorithm
#' @param lambda grids of \eqn{\lambda_1} used in the algorithm
#' @param X_train the training input
#' @param y_train the training output
#' @param X_vali the validation input
#' @param y_vali the validation output
#' @param X_test the test input (if have)
#' @param y_test the test output (if have)
#' @param num_knots num of knots used for the spline basis
#' @param order the order of spline basis
#' @return The sum of \code{x} and \code{y}
#' \item{BB}{the estimated coefficient tensor (in tensor form)}
#' \item{bb}{the estimated intercept term}
#' \item{others}{used for the package developer}
#' @seealso broadcasted_sparsetenreg, sequential_warmstart
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
#' y_train = 1 + crossprod(matrix(BroX_train, c(prod(dim(BB)), n_train)), as.vector(BB)) + signal_level * rnorm(n_train)
#'
#' # the test data
#' X_test = array(runif(prod(c(dim(BB), n_test)), 0, 1), c(dim(BB), n_test))
#' BroX_test = X_test + 0.6 * sin(2 * pi * (X_test - 0.5)^2)
#' y_test = 1 + crossprod(matrix(BroX_test, c(prod(dim(BB)), n_test)), as.vector(BB)) + signal_level * rnorm(n_test)
#'
#'
#' set.seed(1224) # used for validation method
#' n_vali = 0.25 * n_train
#' # generate the test data
#' X_vali = array(runif(prod(c(dim(BB), n_vali)), 0, 1), c(dim(BB), n_test))
#' BroX_vali = X_vali + 0.6 * sin(2 * pi * (X_vali - 0.5)^2)
#' y_vali = 1 + crossprod(matrix(BroX_vali, c(prod(dim(BB)), n_vali)), as.vector(BB)) + signal_level * rnorm(n_vali)
#'
#'
#' R=c(1,2,3,4,5)
#' alpha=c(0,0.5,1)
#' lambda=c(0.01,0.05,0.1,0.5,1,5,10,50,100,500,1000)
#'
#' #res_valiandtest <- validation_broadcasted_sparsetenreg(R,alpha,lambda,X_train,y_train,X_vali,y_vali,X_test,y_test, num_knots=5, order=4)



validation_broadcasted_sparsetenreg<-function(R,alpha,lambda,X_train,y_train,X_vali,y_vali,X_test,y_test, num_knots=5, order=4){

  n=length(y_train)
  #basis transformation
  knots = stats::quantile(c(X_train), probs = c(seq(0, 1, 1/(num_knots - 1))))
  tildePhiX_train = tildePhiX_trans(X_train, knots, order)
  rm(X_train)
  cat=("basis transformation \n")

  tildePhiX_vali=tildePhiX_trans(X_vali, knots, order)
  rm(X_vali)

  tildePhiX_test=tildePhiX_trans(X_test, knots, order)
  rm(X_test)

  cat("basis transformation is completed \n")



  num_R=length(R)
  num_lambda=length(lambda)
  num_alpha=length(alpha)
  num_all=num_R*num_lambda*num_alpha

  b_validation_test_lambda_R_nonlinear=matrix(0,num_all,6)
  BBbig_nonlinear=list()
  betabig <- list()


  for(i1 in 1:num_R){
    for(i2 in 1:num_alpha){
      for(i3 in 1:num_lambda){

        if(i3==1){

          fit=try(quiet(broadcasted_sparsetenreg(TolFun=0.0001,rescale_all=0,rescale_l2=1,r=R[i1],tildePhiX_train,y_train,penalty=c("L1"),knots=NA,restriction = 0,lambda=lambda[i3],gamma=0,alpha_gamma=0,
                                                 alpha=alpha[i2],Replicates = 1,family="gaussian",MaxIter = 10000,manifold=0,warmstart=1, QCQP=1)))
          if("try-error" %in% class(fit)){
          }else{
            fit_nonlinear=fit
          }



          #fit_nonlinear=quiet(broadcasted_sparsetenreg(TolFun=0.0001,rescale_all=0,rescale_l2=1,r=R[i1],tildePhiX_train,y_train,penalty=c("L1"),knots=NA,restriction = 0,lambda=lambda[i3],gamma=0,alpha_gamma=0,
          #                                           alpha=alpha[i2],Replicates = 1,family="gaussian",MaxIter = 10000,manifold=0,warmstart=1, QCQP=1))
        }else{
          fit=try(quiet(broadcasted_sparsetenreg(TolFun=0.0001,rescale_all=0,rescale_l2=1,r=R[i1],tildePhiX_train,y_train,lambda=lambda[i3],gamma=0,alpha_gamma=0,
                                                 alpha=alpha[i2],Replicates = 1,family="gaussian",MaxIter = 10000,manifold=0,warmstart=0, QCQP=1,B0=fit_nonlinear$beta,beta0 =fit_nonlinear$beta0)))

          if("try-error" %in% class(fit)){
          }else{
            fit_nonlinear=fit
          }



          #fit_nonlinear=quiet(broadcasted_sparsetenreg(TolFun=0.0001,rescale_all=0,rescale_l2=1,r=R[i1],tildePhiX_train,y_train,lambda=lambda[i3],gamma=0,alpha_gamma=0,
          #                                         alpha=alpha[i2],Replicates = 1,family="gaussian",MaxIter = 10000,manifold=0,warmstart=0, QCQP=1,B0=fit_nonlinear$beta,beta0 =fit_nonlinear$beta0))

        }

        num_location=(i1-1)*num_alpha*num_lambda+(i2-1)*num_lambda+i3
        cat(stringr::str_c('Complete ','%',round(100*num_location/num_all,2)), "\n")




        b0=fit_nonlinear$beta0
        BB=full_R(fit_nonlinear$beta)
        betabig[[num_location]]=fit_nonlinear$beta

        MSE_nonlinear_validation=prediction_function_nonlinear(b0,BB,tildePhiX_vali,y_vali,family='gaussian')
        #MSE_nonlinear_validation=predict_BroadcasTR(b0,BB,tildePhiX_vali,y_vali,family='gaussian')
        MSE_nonlinear_test=prediction_function_nonlinear(b0,BB,tildePhiX_test,y_test,family='gaussian')
        #MSE_nonlinear_test=predict_BroadcasTR(b0,BB,tildePhiX_test,y_test,family='gaussian')

        BBbig_nonlinear[[num_location]]=BB
        b_validation_test_lambda_R_nonlinear[num_location,1]=b0
        b_validation_test_lambda_R_nonlinear[num_location,2]=MSE_nonlinear_validation
        b_validation_test_lambda_R_nonlinear[num_location,3]=MSE_nonlinear_test
        b_validation_test_lambda_R_nonlinear[num_location,4]=R[i1]
        b_validation_test_lambda_R_nonlinear[num_location,5]=alpha[i2]
        b_validation_test_lambda_R_nonlinear[num_location,6]=lambda[i3]

      }
    }
  }

  #avoid numerical problems
  for(num_location in 1:num_all){
    if(is.nan(b_validation_test_lambda_R_nonlinear[num_location,2])==1){
      b_validation_test_lambda_R_nonlinear[num_location,2]=Inf
    }
    if(is.na(b_validation_test_lambda_R_nonlinear[num_location,2])==1){
      b_validation_test_lambda_R_nonlinear[num_location,2]=Inf
    }
    if(b_validation_test_lambda_R_nonlinear[num_location,2]==0){
      b_validation_test_lambda_R_nonlinear[num_location,2]=Inf
    }

  }

  index_best=which.min(b_validation_test_lambda_R_nonlinear[,2])

  res=list(knots_used = knots, betabig=betabig, BB=BBbig_nonlinear[[index_best]],bb=b_validation_test_lambda_R_nonlinear[index_best,1],MSE_vali=b_validation_test_lambda_R_nonlinear[index_best,2],MSE_pre=b_validation_test_lambda_R_nonlinear[index_best,3],knots=knots,BBbig_nonlinear=BBbig_nonlinear,b_validation_test_lambda_R_nonlinear=b_validation_test_lambda_R_nonlinear)

  return(res)
}

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}
