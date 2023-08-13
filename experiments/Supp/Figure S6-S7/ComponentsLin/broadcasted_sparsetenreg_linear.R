#Final broadcasted nonparametric tensor regression with penalty

broadcasted_sparsetenreg<-function(X_sample,
                                      y,
                                      r=2,
                                      restriction=0,
                                      knots=0,
                                      penalty=c("L1L2"),
                                      lambda=0,
                                      gamma=0,
                                      alpha_gamma=0,
                                      alpha=0,
                                      family="gaussian",
                                      Z_sample=NA,
                                      B0=NA,
                                      B0_list=NA,
                                      BB0=NA,
                                      BB0_list=NA,
                                      Replicates=5,
                                      epsilon=1e-7,
                                      MaxIter=1000,
                                      TolFun=0.0001,
                                      Memorysave=0,
                                      rescale_l2=1,
                                      rescale_all=0,
                                      manifold=0,
                                      warmstart=1,
                                      improvement=1,
                                      stoppingrule_penalty=1,
                                      QCQP=1,
                                   beta0=NA,
                                   beta0_vec=NA,
                                   shrink_factor_number=5,
                                   startsize=3
                                   
)
{
  #used in function parametere
  
  ##knots and restriction
  #knots=..:are used for constraints when using Bspline
  #restriction=1, identifiability restriction when using Bspline
  
  #penalty=c("SCAD"),and lambda=.. is the tunning parameter
  ##if penalty \ne c("SCAD), it is elastic net,
  ###Lambda: tuning parameter for CP component (Coefficient of function)
  ###gamma=0, l2 penalty for CP component (Coefficient of function)
  ###=1,
  ###alpha_gamma: tuning parameter for CP component (Coefficient of function)
  ###alpha=0: l2 penalty for coefficients of basis
  ###=1: l1penalty for the coefficients of basis
  
  ####If lambda=gamma, alpha_gamma=alpha, the data is linear, it is equivalent to tensor linear regresion
  
  #family=“gaussian”: linear regression,
  #Z_sample=NA, other predictor.
  #B0=NA, initial value of all coefficients.

  #Replicates=5,  how many initial point will be chosen
  #epsilon=5*1e-7, the  convergence parameter of glmnet function
  #MaxIter=200, the maximum number of iteration
  #TolFun=0.0001, the parameter  used in stopping rule
  #Memorysave=1, the method that will save memory in the process, however,  it will reduce the speed.
  #Separate_update=1, it seems that no use
  
  
  #warmstar=1, use downsize tensor regression as initial point
  #use ridge regression as initial point
  if(alpha!=0&&alpha!=1){
    elastic_indicator=1
  }else{
    elastic_indicator=0
  }
  
  #this is used to improve some converge problems of glmnet function
  control_iteration=0
  
  #fix the scale of the norm of  coefficient  of basis to be 1. Although the scale dose not matter theoritically, but it have effects nummerically
  Scaleofbasis=1
  
  #this is ued for the warmstart 1, i.e., tensor regression without penalty, but may have downsize
  shrink_factor_scale=5
  
  #this is used for feasible adjustment for QCQP or resacle to get a smaller objective function with penalty
  norm1_iterationstart=1
  
  
  # MaxIter=1000
  #  TolFun=0.0001
  n=length(y)
  p1=dim(X_sample)
  p=p1[-length(p1)]
  
  #used to check if converge
  Convergence_Index=0
  
  
  if(is.na(Z_sample)==1)
  {X=rep(1,n)}else{
    z2=dim(Z_sample)[2]
    X=matrix(0,n,z2+1)
    X[,1:z2]=Z_sample
    X[,z2+1]=rep(n,1)
  }
  X=as.matrix(X)
  x2=dim(X)[2]
  
  d=length(p)
  if(class(X_sample)[1]=="Tensor"){
    TM=X_sample
  }else{
    TM=as.tensor(X_sample)
  }
  
  
  
  #warmstartmethod, unfinished
  #for initial tensor or CP components, if warmstart!=0, then it is the initial tensor of warmstart
  if(warmstart==3){
    K=p[length(p)]
    #this warmstart method is following Hua Zhou's method
   
    if(warmstart==3){
      
      wanmstartdowsize_res=warmstartdowsize(warmstart=warmstart,r=r,X_sample=X_sample,shrink_factor_scale=shrink_factor_scale,B0=B0,B0_list=B0_list,BB0=BB0,BB0_list=BB0_list)
      
      X_sample_reduce=wanmstartdowsize_res$X_sample_reduce
      r_reduce=wanmstartdowsize_res$r_reduce
      B0_reduce=wanmstartdowsize_res$B0_reduce
      B0_list_reduce=wanmstartdowsize_res$B0_list_reduce
      BB0_reduce=wanmstartdowsize_res$BB0_reduce
      BB0_list_reduce=wanmstartdowsize_res$BB0_list_reduce
      shrink_factor=wanmstartdowsize_res$shrink_factor
    
    
    
    BBreduce=broadcasted_tenreg(X_sample=X_sample_reduce, y,r=r_reduce,restriction=0,knots=0,penalty=c("L1L2"),lambda=0,
                                  gamma=0,alpha_gamma=0,alpha=0,family="gaussian",Z_sample=NA,
                                  Replicates=5,epsilon=1e-7,MaxIter=1000,TolFun=0.001,
                                  Memorysave=0,rescale_l2=1,rescale_all=0,
                                  manifold=0,QCQP=0,B0_list = B0_list_reduce,B0=B0_reduce,BB0=BB0_reduce,BB0_list=BB0_list_reduce,beta0=beta0)
    
    testlist=BBreduce
  
    #return(BBreduce)
    beta0=BBreduce$beta0
    BBreduce=full_R(BBreduce$beta)
    BBup=warmstartupsize(warmstart=warmstart,shrink_factor=shrink_factor,p=p,X_sample_reduce_index=X_sample_reduce_index,BBreduce=BBreduce)
    
    
    
    }
    
    warmstartprint(warmstart=warmstart)
    
    cpD=cp(as.tensor(BBup),r_reduce)
    #return(cpD)
    beta=cpD$U
    cp_lambda=cpD$lambdas
    #return(beta)
    for(i in 1:d){
      for(rr in 1:r_reduce){
        betaii=beta[[i]]  
        beta[[i]][,rr]= betaii[,rr]*(cp_lambda[rr])^(1/d)
      }
    }
    
    
    
    
  }else{
    #the new warmstart 1 is the sequential warmtsart.
    if(warmstart==1){
      fit_nonlinear_sequantial_warmstart=sequential_warmstart_study(TolFun=0.01,rescale_all=0,rescale_l2=0,r=r,X_sample=X_sample,y=y,penalty=c("L1"),knots=knots,restriction = 0,lambda=0.00,gamma=0.000,alpha_gamma=0,
                                                                    alpha=0,Replicates = NA,family="gaussian",MaxIter = 50,manifold=0,warmstart=3, QCQP=0,shrink_factor_number=shrink_factor_number,startsize=startsize)
      beta=fit_nonlinear_sequantial_warmstart$bbeta
      beta0=fit_nonlinear_sequantial_warmstart$bbeta0
    }else{
    beta=list()
    }
  }
  
  
  
  
  
  
  
  
  rm(X_sample)
  gc()
  
  Md=list()
  #beta=list()
  #dimTM=dim(TM)
  if(Memorysave==1){
  }else{
    for (dd in 1:d){
      Md[[dd]]=unfold(TM,c(d+1,dd),c(1:d)[-dd])
    }
    rm(TM)
    gc()}
  #use list in R to save some results
  dev_final=Inf
  
  #This  is an updating of hua zhou's code
  dev_final_penalty=Inf
  
  
  #Objective function final
  #up to now, only uesed of l2 penalty
  ob_f_final=0
  ob_f_loss_penalty_final=0
  
  
  #f means without penalty, just loss
  ob_f_differentinitial=c()

  #loss for various initial points(may be used in the future)
  ob_f_loss_penalty_differentinitial=c()
  
  #linearpart part for different initial values
  linearpart_differentinitial=list()
  
  
  
  ##This is  used for debug
  betaold=0
  
  
  #The iteration number 
  num_final=0
  
  #The replicate. If warmstart=0, then it just equals to the lengths of B0_list or BB0_list
  if(warmstart==0){
  if(is.na(sum(B0_list[[1]][[1]]))!=1){
    Replicates=length(B0_list) 
  }
  if(is.na(sum(BB0_list[[1]][[1]]))!=1){
    Replicates=length(BB0_list) 
  }
  }
  
  #If only on initial point, then, the replicate equals to the use's setting. If many, then the replicate equals to many.
  
  
  for (rep in 1:Replicates){
    
    #Objective function final
    
    ob_f=c()
    ob_f_loss_penalty=c()
    #initialize tensor regression coefficients from uniform [-1,1]  or standard norm distribution 
    set.seed(rep*10)
    
    #when use warmstar, we also can use random initial point and compare the final objective function value, and choose the best
    if(is.na(sum(B0_list[[1]][[1]]))!=1||is.na(sum(BB0_list[[1]][[1]]))!=1){
    }else{
    if(rep>1){
    B0=NA
    BB0=NA
    warmstart=0
    }
    }
    
    
    if(is.na(sum(B0[[1]]))==1&is.na(sum(BB0[[1]]))==1&is.na(sum(B0_list[[1]][[1]]))==1&is.na(sum(BB0_list[[1]][[1]]))==1){
      if(warmstart!=0){
        if(QCQP==1||manifold==1){
        #to fix the  norm 1
        for(i in 1:r){
          norm1scale=sqrt(sum(beta[[d]][,i]^2))*Scaleofbasis
        beta[[d]][,i]=beta[[d]][,i]/norm1scale
        beta[[1]][,i]=beta[[1]][,i]*norm1scale
        }
        }
        
      }else{
        for(i in 1:d){
          #both partern is OK
          beta[[i]]=matrix(runif(p[i]*r,-1,1),p[i],r)
          #beta[[i]]=matrix(rnorm(p[i]*r),p[i],r)
        }
        if(QCQP==1||manifold==1){
          #to fix the  norm 1
          for(i in 1:r){
            norm1scale=sqrt(sum(beta[[d]][,i]^2))*Scaleofbasis
            beta[[d]][,i]=beta[[d]][,i]/norm1scale
            beta[[1]][,i]=beta[[1]][,i]*norm1scale
          }
        }
      }
    }else{
      if(warmstart!=0){
      #If there is a warmstart, then just use the  warmstart initial point  
      }else{
      if(is.na(sum(B0[[1]]))!=1){
      for(i in 1:d){
        beta[[i]]=matrix(B0[[i]],p[i],r)
      }
      }
      if(is.na(sum(BB0[[1]]))!=1){
        print("CP decomposition for the initial tensor")
        cpBB0=cp(as.tensor(BB0),r)
        #return(cpD)
        beta=cpBB0$U
        cp_lambda=cpBB0$lambdas
        #return(beta)
        for(i in 1:d){
          for(rr in 1:r){
            betaii=beta[[i]]  
            beta[[i]][,rr]= betaii[,rr]*(cp_lambda[rr])^(1/d)
          }
        }
        
      }
      
      if(is.na(sum(B0_list[[1]][[1]]))!=1){
        beta=B0_list[[rep]]
      }
      if(is.na(sum(BB0_list[[1]][[1]]))!=1){
        BB0_middle=BB0_list[[rep]]
        print("CP decomposition for the initial tensor")
        cpBB0=cp(as.tensor(BB0_middle),r)
        #return(cpD)
        beta=cpBB0$U
        cp_lambda=cpBB0$lambdas
        #return(beta)
        for(i in 1:d){
          for(rr in 1:r){
            betaii=beta[[i]]  
            beta[[i]][,rr]= betaii[,rr]*(cp_lambda[rr])^(1/d)
          }
        }
        
      }
      }
      
      if(QCQP==1||manifold==1){
        #to fix the  norm 1
        for(i in 1:r){
          norm1scale=sqrt(sum(beta[[d]][,i]^2))*Scaleofbasis
          beta[[d]][,i]=beta[[d]][,i]/norm1scale
          beta[[1]][,i]=beta[[1]][,i]*norm1scale
        }
      }
      
      #Replicates=1
    }  
    


    
    #main loop 
    for (iter in 1:MaxIter){
      num=iter
      if(iter==1){
        #actually, this is the initial value of u
        if(is.na(sum(beta0))==1){
        if(warmstart!=0){
        }else{
        df<-data.frame(X,y)
        fit<-glm(y~.-1,family = family,data=df,control = list(maxit = 1000000,epsilon = 1e-8))
        beta0=coef(fit)
        }
        }else{
            beta0=beta0
          }
        if(is.na(sum(beta0_vec))!=1){
          beta0=beta0_vec[rep]
        }
        
        
        linearpart=as.matrix(X)%*%beta0
        dev0=loglike(y,linearpart,family=family)
        
        #the objective function and penalized objective function
        #only uesd when l2 penalty
        
        #this is an updating 
        
        
        
        dev0_penalty=dev0+n*lambda*penaltyvalue(beta=beta,alpha=alpha,manifold=manifold)
        
        
        
        #ob_f[1]=dev0
        #ob_f_loss_penalty[1]=ob_f[1]+dev0_penalty
#        if(is.na(sum(B0[[1]]))!=1){
#          ob_f[1]=dev0
#          ob_f_loss_penalty[1]=ob_f[1]+dev0_penalty
            
#        }else{
        ob_f[1]=0
        ob_f_loss_penalty[1]=0
#        }
      }else{
        #Stopping rule is in this part, which means  iteration >1
        eta=Xj%*%as.vector(beta[[d]])
        
        X_eta=matrix(0,n,x2+1)
        X_eta[,1:x2]=X
        X_eta[,x2+1]=eta
        
        df<-data.frame(X_eta,y)
        #Forget the mean of Separate_update.
        #Although, I forget this, i can use offset, to improve hua zhou's code
        {
          #used to avoid the numerical problems
          if(sum(abs(eta))<1e-16){
            print("lambda is too big, the coefficients of CP component close to 0")
            print(lambda)
            print(r)
            if(family=='binomial'){
              pre_error=sum(abs(as.numeric(linearpart>0)-y))/n
            }else{
              #pre_error is distance here
              pre_error=norm(as.matrix(y)-linearpart,'f')^2/n
            }
            
            
            
            
            
            #The sparse parameters
            pzeros=0
            for(i in 1:d){
                pzeros=pzeros+sum(beta[[i]]==0)
            }
            
            
            
            #output BIC of the final model. Note deviance = -2*log-likelihood
            if(d==2){
                df=r*(p[1]+p[2]-r)+x2-pzeros
                #   BIC=log(dev_final)+log(n)*(r*(p[1]+p[2]-r)+x2)
            }else{
                df=r*(sum(p[1:(d-1)])-d+1)+x2-pzeros
                #  BIC=log(dev_final)+log(n)*(r*(sum(p[1:(d-1)])-d+1)+x2);
            }
            BIC=n*(log(2*pi)+log(dev0))+n+log(n)*df
            
            tunning_index=list(deviance=n*(log(2*pi)+log(dev0))+n,df=df)
            
            early_stop_for_lambda=1
            

            
            
            test=list(beta0=beta0,beta=beta,dev=dev0,error=pre_error,num=num,ob_f_final=ob_f/n,ob_f_loss_penalty_final=ob_f_loss_penalty/n,BIC=BIC,tunning_index=tunning_index,early_stop_for_lambda=early_stop_for_lambda)

            return(test)
            
            
            
            
            
            
            
            
            
            #        test=list(beta0=beta0,beta=beta,dev=dev0,error=pre_error,num=num)
            #return(test)
          }
          
          if(improvement==1&&iter> control_iteration){
            #print('haha')
            df<-data.frame(X,y)
            #use  offset or just  minus it in y?
            fit<-glm(y~.-1,family = family,offset=eta,data=df,control = list(maxit = 10000000,epsilon = epsilon))
            beta_eta=c(coef(fit),1)
            end=length(beta_eta)
            beta0=beta_eta[-end] 
            
          }else{
            #print('hahahahahahahahahahhahahahahhahaha')
            fit<-glm(y~.-1,family = family,data=df,control = list(maxit = 10000000,epsilon = epsilon))
            beta_eta=coef(fit)
            end=length(beta_eta)
            beta0=beta_eta[-end]
          }
        }
        
        #Zy: Xj also need to be rescale test
        ########
        
        
        ########
        
        
        
        
        
        ####this below is hua zhou's code position, definitely, this is the loss function,
        #stopping rule
        if(improvement==1&&iter> control_iteration){
          linearpart=eta+beta0
          devtmp=loglike(y,linearpart,family=family) 
          
        }else{
          linearpart=X_eta%*%beta_eta
          devtmp=loglike(y,linearpart,family=family) 
        }
        

        
        
        ################an update to hua zhou's stopping rule
        if(improvement==1&&iter> control_iteration){
          #then updating the intercept, do not update the scale of CP components.
          

          if(iter>norm1_iterationstart){
          #feasible adjustment and fix norm 1
          #if it is infeasible, then it will be feasible
          #if it if feasible,the the objective funtion (with penalty) will be smaller
          if(QCQP==1||manifold==1){  
#          if(norm(beta[[d]],'f')^2>r){
            #to fix the  norm 1
            for(i in 1:r){
              norm1scale=sqrt(sum(beta[[d]][,i]^2))*sqrt(Scaleofbasis)
              if(norm1scale>1e-2){
              beta[[d]][,i]=beta[[d]][,i]/norm1scale
              beta[[1]][,i]=beta[[1]][,i]*norm1scale
              }
            }
#          }
          }
          }
          
#          before=penaltyvalue(beta=beta,alpha=alpha,manifold=manifold,QCQP=QCQP)
          beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)  
#          after=penaltyvalue(beta=beta,alpha=alpha,manifold=manifold,QCQP=QCQP)
#          print(before-after)
        }else{
          beta=arrange(beta,beta_eta[length(beta_eta)]) 
          
#          before=penaltyvalue(beta=beta,alpha=alpha,manifold=manifold,QCQP=QCQP)
          beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)
          
#          after=penaltyvalue(beta=beta,alpha=alpha,manifold=manifold,QCQP=QCQP)
#          print(before-after)
        }
        #########################################
        
        
        
        #the objective values
        ob_f[(d+1)*(iter-1)+1]=devtmp
        ob_f_loss_penalty[(d+1)*(iter-1)+1]=devtmp+n*lambda*penaltyvalue(beta=beta,alpha=alpha,manifold=manifold,QCQP=QCQP)
        ####this above is hua zhou's code position
        
        #This is an updating of stopping rule
        devtmp_penalty=devtmp+n*lambda*penaltyvalue(beta=beta,alpha=alpha,manifold=manifold,QCQP=QCQP)
        
        #print(devtmp_penalty)
        
        

        
        ###############used for boundary beyound, devtmp, linearpart too big or, linearpart includes NA?
        if(is.nan(devtmp)==1 || is.na(devtmp)==1||is.na(devtmp_penalty)==1||is.nan(devtmp_penalty)==1){
          #In the former iteration, it has been reach the boundary
          print("error type 1")
          print('beyond the number boundary in the software, devtmp not a number, or, linearpart NA')
          print(iter)
          print(r)
          print(lambda)
          if(family=='binomial'){
            pre_error=sum(abs(as.numeric(linearpart>0)-y))/n
          }else{
            #pre_error is distance here
            pre_error=norm(as.matrix(y)-linearpart,'f')^2/n
          }
          
          
          #The sparse parameters
          pzeros=0
          for(i in 1:d){
              pzeros=pzeros+sum(beta[[i]]==0)
          }
          
          
          
          #output BIC of the final model. Note deviance = -2*log-likelihood
          if(d==2){
              df=r*(p[1]+p[2]-r)+x2-pzeros
              #   BIC=log(dev_final)+log(n)*(r*(p[1]+p[2]-r)+x2)
          }else{
              df=r*(sum(p[1:(d-1)])-d+1)+x2-pzeros
              #  BIC=log(dev_final)+log(n)*(r*(sum(p[1:(d-1)])-d+1)+x2);
          }
          BIC=n*(log(2*pi)+log(dev0))+n+log(n)*df
          
          tunning_index=list(deviance=n*(log(2*pi)+log(dev0))+n,df=df)
          
          
          
          
          
          
          
          
           early_stop_for_lambda=1
          
          test=list(beta0=beta0,beta=beta,dev=dev0,error=pre_error,num=num,ob_f_final=ob_f/n,ob_f_loss_penalty_final=ob_f_loss_penalty/n,BIC=BIC,tunning_index=tunning_index, early_stop_for_lambda=early_stop_for_lambda)
          return(test)   
        }
        
        
        ##############
        diffdev=devtmp-dev0
        dev0=devtmp
        
        
        ##This is an updating to hua zhou's stopping rule
        diffdev_penalty=devtmp_penalty-dev0_penalty
        dev0_penalty=devtmp_penalty
        if(stoppingrule_penalty==1)
        {
          v1=abs(diffdev_penalty)
          #print(v1)
          v2=TolFun*(abs(dev0_penalty)+1)
          #print(v2)
        }else{
          v1=abs(diffdev)
          #print(v1)
          v2=TolFun*(abs(dev0)+1)
          #print(v2)
        }
        
    
        #This is an updating to hua zhou's code

        if(v1<v2){
          #print(v1)
          #print(v2)
          Convergence_Index=1
          #fix==1 means fit the some entries to be one, refer to wiki of raymond
         {
            #This is a improvement of hua zhou's code
            #when the  iteration number big, just update the intercept, 
            #when iteration number is not big, for numerical reason, update them togeter
            if(improvement==1&&iter> control_iteration){
              beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)  
            }else{
              beta=arrange(beta,beta_eta[length(beta_eta)]) 
              beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)
            }
          }
          print("reach the stopping rule")
          break
        }
        ##ZY thinks this following should be placed in above.
        #update scale of array coefficients and standardize(Old version)(339-359)
        {
          #this is  an improvement to hua zhou's code
          if(improvement==1&&iter> control_iteration){
            beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)
          }else{
            beta=arrange(beta,beta_eta[length(beta_eta)])
            beta=rescale(beta,rescale_all = rescale_all,rescale_l2 = rescale_l2,elastic=elastic_indicator,alpha=alpha)
          }
        }
      }
      
      # cyclic update of the array coefficients
      eta0=as.matrix(X)%*%beta0
      #loop for CP components
      for(j in 1:d){
        if (j==1){
          cumkr=t(as.matrix(rep(1,r)))
        }
        if (j==d){
          #same memory method
          if(Memorysave==1){
            Md[[j]]=unfold(TM,c(d+1,j),c(1:d)[-j]) 
            Xj=matrix(Md[[j]]@data%*%cumkr,n,p[j]*r)
            Md[[j]]=0
            gc()
          }else{
            ##this  is used for debug
            #            Xj_old=Xj
            ####
            Xj_old=Xj
            #            
            Xj=matrix(Md[[j]]@data%*%cumkr,n,p[j]*r)
          }
        }else{
          if(j==(d-1)){
            if(Memorysave==1){
              Md[[j]]=unfold(TM,c(d+1,j),c(1:d)[-j]) 
              Xj=matrix(Md[[j]]@data%*%khatri_rao(beta[[d]],cumkr),n,p[j]*r)
              Md[[j]]=0
              gc()
            }else{
              ##this  is used for debug
              #             Xj_old=Xj
              ####
              Xj=matrix(Md[[j]]@data%*%khatri_rao(beta[[d]],cumkr),n,p[j]*r)
            }
          }else{
            if(Memorysave==1){
              Md[[j]]=unfold(TM,c(d+1,j),c(1:d)[-j]) 
              Xj=matrix(Md[[j]]@data%*%khatri_rao(khatri_rao_list(beta[c(d:(j+1))]),cumkr),n,p[j]*r)
              Md[[j]]=0
              gc()
            }else{
              ##this  is used for debug
              #              Xj_old=Xj
              ####
              Xj=matrix(Md[[j]]@data%*%khatri_rao(khatri_rao_list(beta[c(d:(j+1))]),cumkr),n,p[j]*r)}
          }
        }
        
        Xj_eta=matrix(0,n,p[j]*r+1)
        #Corresponding to the CP components
        Xj_eta[,1:(p[j]*r)]=Xj
        #Corresponding to the intercept
        Xj_eta[,p[j]*r+1]=eta0
        #print(Xj_eta)
        #########################################used for boundary beyound, Xj
        if(sum(is.nan(Xj))!=0){
          #In the former iteration, it has been reach the boundary
          print("error type 2")
          print('beyond the number boundary in the software, Xj includes NAN')
          print(iter)
          print(r)
          print(lambda)
          if(family=='binomial'){
            pre_error=sum(abs(as.numeric(linearpart>0)-y))/n
          }else{
            #pre_error is distance here
            pre_error=norm(as.matrix(y)-linearpart,'f')^2/n
          }
          
          
          #The sparse parameters
          pzeros=0
          for(i in 1:d){
              pzeros=pzeros+sum(beta[[i]]==0)
          }
          
          
          
          #output BIC of the final model. Note deviance = -2*log-likelihood
          if(d==2){
              df=r*(p[1]+p[2]-r)+x2-pzeros
              #   BIC=log(dev_final)+log(n)*(r*(p[1]+p[2]-r)+x2)
          }else{
              df=r*(sum(p[1:(d-1)])-d+1)+x2-pzeros
              #  BIC=log(dev_final)+log(n)*(r*(sum(p[1:(d-1)])-d+1)+x2);
          }
          BIC=n*(log(2*pi)+log(dev0))+n+log(n)*df
          
          tunning_index=list(deviance=n*(log(2*pi)+log(dev0))+n,df=df)
          
          
          
          early_stop_for_lambda=1
          
          test=list(beta0=beta0,beta=beta,dev=dev0,error=pre_error,num=num,ob_f_final=ob_f/n,ob_f_loss_penalty_final=ob_f_loss_penalty/n,BIC=BIC,tunning_index=tunning_index,early_stop_for_lambda=early_stop_for_lambda)
          return(test)  
        }
        
        if(sum(Xj)==0){
          #In the former iteration, it has been reach the boundary
          print("error type 3")
          print('lambda is too big, Xj (corresponding to the coefficients of CP components) is 0')
          print(iter)
          print(r)
          print(lambda)
          if(family=='binomial'){
            pre_error=sum(abs(as.numeric(linearpart>0)-y))/n
          }else{
            #pre_error is distance here
            pre_error=norm(as.matrix(y)-linearpart,'f')^2/n
          }
          
          
          #The sparse parameters
          pzeros=0
          for(i in 1:d){
              pzeros=pzeros+sum(beta[[i]]==0)
          }
          
          
          
          #output BIC of the final model. Note deviance = -2*log-likelihood
          if(d==2){
              df=r*(p[1]+p[2]-r)+x2-pzeros
              #   BIC=log(dev_final)+log(n)*(r*(p[1]+p[2]-r)+x2)
          }else{
              df=r*(sum(p[1:(d-1)])-d+1)+x2-pzeros
              #  BIC=log(dev_final)+log(n)*(r*(sum(p[1:(d-1)])-d+1)+x2);
          }
          BIC=n*(log(2*pi)+log(dev0))+n+log(n)*df
          
          tunning_index=list(deviance=n*(log(2*pi)+log(dev0))+n,df=df)
          
          
          
          early_stop_for_lambda=1
          
          test=list(beta0=beta0,beta=beta,dev=dev0,error=pre_error,num=num,ob_f_final=ob_f/n,ob_f_loss_penalty_final=ob_f_loss_penalty/n,BIC=BIC,tunning_index=tunning_index,early_stop_for_lambda=early_stop_for_lambda)
# to avoid that some cases.
#    break;
          return(test)
        }
        #########################################
        
        
        
        
        
        
        
          if(j<d){
            #for SCAD penalty
            if(penalty==c("SCAD")){
              df_beta=data.frame(Xj_eta=Xj_eta,y=y)
              fit<-glmreg(y~.,standardize=FALSE,data=df_beta,family=c("gaussian"),penalty=c("enet"),lambda=lambda)
              Bd_new_vec_and=coef(fit,exact=TRUE)[-1]}
            
            #for L1 or L2 penalty
            if(penalty!=c("SCAD")){
                
                if(improvement==1&&iter>control_iteration) 
                {
                  #print(j)
                  fit<-glmnet(Xj_eta[,-(p[j]*r+1)],y-Xj_eta[,(p[j]*r+1)],family=family,lambda=lambda,intercept=1,alpha = alpha,maxit=100000000,standardize = FALSE,thresh=epsilon)
                  Bd_new_vec_and=coef(fit,exact=TRUE)[-1] 
                  ##This is used for debug
                  betaold=beta
                  ##
                 
                  
                  beta[[j]]=matrix(Bd_new_vec_and,p[j],r)
                  eta0=Xj_eta[,(p[j]*r+1)]+coef(fit,exact=TRUE)[1]
                  
                  linearpart=Xj_eta[,-(p[j]*r+1)]%*%Bd_new_vec_and+eta0
                  devtmp=loglike(y,linearpart,family=family)
                  #print(devtmp)
 
                }else{
                  #This is the current function
                  ###a test
                  p.fac = rep(1, p[j]+1)
                  p.fac[p[j]+1]=0
                  
                  #intercept is in (p[j]+1)-th colume
                  fit<-glmnet(Xj_eta,y,family=family,lambda=lambda,intercept=0,penalty.factor = p.fac,alpha = alpha,maxit=100000000,standardize = FALSE,thresh=epsilon)
                  Bd_new_vec_and=coef(fit,exact=TRUE)[-1] 
                  
                  #########################################used for setting checking,  Bj j<d
                  if(sum(abs(Bd_new_vec_and[-length(Bd_new_vec_and)]))==0){
                    print("error type 5")
                    #In the former iteration, it has been reach the boundary
                    print('beyond the number boundary in the software or lambda is too big, Bj, j<d, return 0')
                    
                    print(iter)
                    print(r)
                    print(lambda)
                    if(family=='binomial'){
                      pre_error=sum(abs(as.numeric(linearpart>0)-y))/n
                    }else{
                      #pre_error is distance here
                      pre_error=norm(as.matrix(y)-linearpart,'f')^2/n
                    }
                    
                    #The sparse parameters
                    pzeros=0
                    for(i in 1:d){
                        pzeros=pzeros+sum(beta[[i]]==0)
                    }
                    
                    
                    
                    #output BIC of the final model. Note deviance = -2*log-likelihood
                    if(d==2){
                        df=r*(p[1]+p[2]-r)+x2-pzeros
                        #   BIC=log(dev_final)+log(n)*(r*(p[1]+p[2]-r)+x2)
                    }else{
                        df=r*(sum(p[1:(d-1)])-d+1)+x2-pzeros
                        #  BIC=log(dev_final)+log(n)*(r*(sum(p[1:(d-1)])-d+1)+x2);
                    }
                    BIC=n*(log(2*pi)+log(dev0))+n+log(n)*df
                    
                    tunning_index=list(deviance=n*(log(2*pi)+log(dev0))+n,df=df)
                    
                    
                    
                    early_stop_for_lambda=1
                    
                    
                    test=list(beta0=beta0,beta=beta,dev=dev0,error=pre_error,ob_f_final=ob_f/n,ob_f_loss_penalty_final=ob_f_loss_penalty/n,BIC=BIC,tunning_index=tunning_index,early_stop_for_lambda=early_stop_for_lambda)
                    return(test)  
                  }
                  ###############################################################################
                  
                  beta[[j]]=matrix(Bd_new_vec_and[-length(Bd_new_vec_and)],p[j],r)
                  eta0=eta0*Bd_new_vec_and[length(Bd_new_vec_and)]  
                  
                  #The objective values
                  linearpart=Xj_eta%*%Bd_new_vec_and+eta0
                  devtmp=loglike(y,linearpart,family=family)
                 
                  
                }
              }
            

            ob_f[(d+1)*(iter-1)+1+j]=devtmp
            devtmp_penalty=devtmp+n*lambda*penaltyvalue(beta=beta,alpha=alpha,manifold=manifold,QCQP=QCQP)
            ob_f_loss_penalty[(d+1)*(iter-1)+1+j]=devtmp_penalty
            if(iter>1){
              if(ob_f_loss_penalty[(d+1)*(iter-1)+1+j]>ob_f_loss_penalty[(d+1)*(iter-1)+j]){
               # print(c(iter,j))
              }
            }
            
            
          }else{
            if(QCQP==1||manifold==1){
            if(QCQP==1){
              #there are some other methods that  works not good for the QCQP
              #beta[[j]]=coodinateforalpha(y=y-Xj_eta[,(p[j]*r+1)],DX=Xj_eta[,-(p[j]*r+1)],alphaold=beta[[j]],alphaTolFun=0.00001,Maxiterforalpha=2000)
              #fit=Socpforalpha(y=y-Xj_eta[,(p[j]*r+1)],DX=Xj_eta[,-(p[j]*r+1)],alphaold=beta[[j]])
              #beta[[j]]=matrix(fit[-(p[j]*r+1)],p[j],r)
              
              #Scaleofbasis=1
              
              betajstandard=beta[[j]]
              linearpartjstandard=Xj_eta[,-(p[j]*r+1)]%*%c(beta[[j]])+Xj_eta[1,(p[j]*r+1)]
              devtmpstandard=loglike(y,linearpartjstandard,family=family)
              
              
              fit=dualascentforalpha(y=y-Xj_eta[,(p[j]*r+1)],DX=Xj_eta[,-(p[j]*r+1)],alphaold=beta[[j]],C=Scaleofbasis,Maxiterfordual=1000,TolFundualascent=0.001,MaxIterInner=40)
              beta[[j]]=matrix(fit$alphavec_final,p[j],r)
             #  hehe=Socpforalpha(y=y-Xj_eta[,(p[j]*r+1)],DX=Xj_eta[,-(p[j]*r+1)],alphaold=beta[[j]],C=1)
             #  beta[[j]]=matrix(hehe[-length(hehe)],p[j],r)
 
              linearpart=Xj_eta[,-(p[j]*r+1)]%*%c(beta[[j]])+Xj_eta[1,(p[j]*r+1)]
              devtmp=loglike(y,linearpart,family=family)
              
              if(devtmp>devtmpstandard){
                beta[[j]]=betajstandard
                linearpart=linearpartjstandard
                devtmp=devtmpstandard
              }
              
              
              
              ob_f[(d+1)*(iter-1)+1+j]=devtmp
              devtmp_penalty=devtmp+n*lambda*penaltyvalue(beta,alpha,manifold=manifold,QCQP=QCQP)
              ob_f_loss_penalty[(d+1)*(iter-1)+1+j]=devtmp_penalty
              
              
              
            }else{
                beta_test=beta[[j]]
                betaold=beta[[j]]
                
                eta0=Xj_eta[,(p[j]*r+1)]

                #print(j)
                if(iter==1){
                  Manifoldresult=update_alpha(Y=y-eta0,BD1=betaold,D=Xj_eta[,-(p[j]*r+1)],K=p[j],r=r,maxupdate_out=100,maxupdate_in=60,scalefactor=0.5,thresholdformani=1e-9)
             
                  beta[[j]]=Manifoldresult$Xl
                  alphamanifold=Manifoldresult$alphafast
                  fl=Manifoldresult$fl
                }else{
                  Manifoldresult=update_alpha(fast=0,alphafast=alphamanifold,Y=y-eta0,BD1=betaold,D=Xj_eta[,-(p[j]*r+1)],K=p[j],r=r,maxupdate_out=1,maxupdate_in=60,scalefactor=0.5,thresholdformani=1e-9)
                  beta[[j]]=Manifoldresult$Xl
                  alphamanifold=Manifoldresult$alphafast
                  fl=Manifoldresult$fl
                }
                
                #The objective values
                linearpart=Xj_eta[,-(p[j]*r+1)]%*%(c(beta[[j]]))+eta0
                devtmp=loglike(y,linearpart,family=family)

                #the objective values
                ob_f[(d+1)*(iter-1)+1+j]=devtmp
                devtmp_penalty=devtmp+n*lambda*penaltyvalue(beta,alpha,manifold=manifold,QCQP=QCQP)
                ob_f_loss_penalty[(d+1)*(iter-1)+1+j]=devtmp_penalty
                
                
                
              }
              }else{
            if(restriction==1){
                  lbound=knots[1]
                  rbound=knots[length(knots)]
                  interval_knots=knots[2:(length(knots)-1)]
                  uu=ibs(rbound, knots = interval_knots, degree = 3, intercept = TRUE,Boundary.knots = c(lbound,rbound))
                  con=matrix(0,r,r*p[[d]])
                  for (rr in 1:r){
                    con[rr,((rr-1)*p[[d]]+1):(rr*p[[d]])] =uu 
                  }
                  #return(matrix(Xj_eta[,1:p[[d]]*r],n,p[[d]]*r))
                  X_aug=rbind(matrix(Xj_eta[,1:p[[d]]*r],n,p[[d]]*r),alpha_gamma*diag(p[[d]]*r))
                  XTX=t(X_aug)%*%X_aug
                  YTX=t(as.matrix(c(y- Xj_eta[,p[d]*r+1],rep(0,p[[d]]*r))))%*%X_aug
                  #heheda=list(XTX=XTX,con=con,YTX=YTX)
                  #return(heheda)
                  fit=solve.QP(Dmat=XTX,dvec=YTX,Amat=t(con),bvec=rep(0,r),meq = r)
                  
                  Bd_new_vec_and=c(fit$solution, Xj_eta[1,p[d]*r+1])
                  
                  
                  #also  update the intercept, scaling  and addition.
                  #dfres=data.frame(Xj_eta=Xj_eta,y=y)
                  ##print(sum(is.na(Xj_eta)))
                  ##print(y)
                  ##print(con)
                  #modelres=lm(y~.-1,data=dfres)
                  ##return(dfres)
                  ##print(modelres)
                  #fit<-conLM.lm(object = modelres, constraints = con, rhs = rep(0,r), neq = r)
                  #Bd_new_vec_and=coef(fit,exact=TRUE)
                }else{
                  
                  if(improvement==1&&iter>control_iteration) 
                  {
                    #print(j)
                    fit<-glmnet(Xj_eta[,-(p[j]*r+1)],y-Xj_eta[,(p[j]*r+1)],family=family,lambda=gamma,intercept=1,alpha = alpha_gamma,maxit=10000000,standardize = FALSE,thresh=epsilon)
                    Bd_new_vec_and=coef(fit,exact=TRUE)[-1] 
                    ##This is used for debug
                    #betaold=beta
                    ##
                    
                    beta[[j]]=matrix(Bd_new_vec_and,p[j],r)
                    eta0=Xj_eta[,(p[j]*r+1)]
                    
                    linearpart=Xj_eta[,-(p[j]*r+1)]%*%Bd_new_vec_and+eta0+coef(fit,exact=TRUE)[1] 
                    devtmp=loglike(y,linearpart,family=family)
                    
                  }else{
                    p.fac = rep(1, p[j]+1)
                    p.fac[p[j]+1]=0
                    
                    fit<-glmnet(Xj_eta,y,family=family,lambda=gamma,intercept=0,penalty.factor = p.fac,alpha = alpha_gamma,maxit=10000000,standardize = FALSE,thresh=epsilon)
                    Bd_new_vec_and=coef(fit,exact=TRUE)[-1]
                    
                    
                    beta[[j]]=matrix(Bd_new_vec_and[-length(Bd_new_vec_and)],p[j],r)
                    ## eta0=eta0*Bd_new_vec_and[length(Bd_new_vec_and)]+coef(fit,exact=TRUE)[1]
                    eta0=eta0*Bd_new_vec_and[length(Bd_new_vec_and)]
                    linearpart=Xj_eta%*%Bd_new_vec_and+eta0
                    devtmp=loglike(y,linearpart,family=family)
                    #              print(j)
                    #             print(devtmp)
                  }
                  
                }
              }
                
                
                #this is  used for "double" penalized
                ########################################used for boundary beyound, Bj, j=d
                if(sum(abs(Bd_new_vec_and[-length(Bd_new_vec_and)]))==0){
                  print("error type 6")
                  #In the former iteration, it has been reach the boundary
                  print('beyond the number boundary in the software or gamma is too large, Bd return 0, fix=0')
                  print(iter)
                  print(r)
                  print(lambda)
                  if(family=='binomial'){
                    pre_error=sum(abs(as.numeric(linearpart>0)-y))/n
                  }else{
                    #pre_error is distance here
                    pre_error=norm(as.matrix(y)-linearpart,'f')^2/n
                  }
                  
                  
                  #The sparse parameters
                  pzeros=0
                  for(i in 1:d){
                      pzeros=pzeros+sum(beta[[i]]==0)
                  }
                  
                  
                  
                  #output BIC of the final model. Note deviance = -2*log-likelihood
                  if(d==2){
                      df=r*(p[1]+p[2]-r)+x2-pzeros
                      #   BIC=log(dev_final)+log(n)*(r*(p[1]+p[2]-r)+x2)
                  }else{
                      df=r*(sum(p[1:(d-1)])-d+1)+x2-pzeros
                      #  BIC=log(dev_final)+log(n)*(r*(sum(p[1:(d-1)])-d+1)+x2);
                  }
                  BIC=n*(log(2*pi)+log(dev0))+n+log(n)*df
                  
                  tunning_index=list(deviance=n*(log(2*pi)+log(dev0))+n,df=df)
                  
                  early_stop_for_lambda=1
                  
                  
                  test=list(beta0=beta0,beta=beta,dev=dev0,error=pre_error,ob_f_final=ob_f/n,ob_f_loss_penalty_final=ob_f_loss_penalty/n,BIC=BIC, tunning_index= tunning_index,early_stop_for_lambda=early_stop_for_lambda)
                  return(test)  
                }
                ################################################################

                
                
                
                #the objective values
                ob_f[(d+1)*(iter-1)+1+j]=devtmp
                devtmp_penalty=devtmp+n*lambda*penaltyvalue(beta,alpha,manifold=manifold,QCQP=QCQP)
                ob_f_loss_penalty[(d+1)*(iter-1)+1+j]=devtmp_penalty

                
              
            }
          
          
        cumkr=khatri_rao(beta[[j]],cumkr)
        
        
      }
      
      
    }  
    
    #the final linearpart for different initial points
    linearpart_differentinitial[[rep]]=linearpart
    
    #the final loss for different initial points.
    ob_f_differentinitial[rep]=dev0
    ob_f_loss_penalty_differentinitial[rep]=dev0_penalty
    
    #This is used  for  the stopping rule
    if(stoppingrule_penalty==1)
    {
      #the loss+penalty 
      v3=dev0_penalty
      v4=dev_final_penalty
    }else{
      #just the loss, used by hua zhou's code
      v3=dev0
      v4=dev_final
    }
    
    

    #this is an updating of hua zhou's code
    # record if it has a smaller deviance
    if(v3<v4)
    {
      beta0_final=beta0
      beta_final=beta
      dev_final=dev0
      linearpart_final=linearpart
      dev_final_penalty=ob_f_loss_penalty[length(ob_f_loss_penalty)]
      #finaliteration number for a specific initial point.
      num_final=num
      
      #get the final objective funtion that we used
      #rescale by the sample size n
      ob_f_final=ob_f/n
      ob_f_loss_penalty_final=ob_f_loss_penalty/n
 
    }
 
  }
  

  #The sparse parameters
  pzeros=0
  for(i in 1:d){
      pzeros=pzeros+sum(beta_final[[i]]==0)
  }
  
  
  
  #output BIC of the final model. Note deviance = -2*log-likelihood
  if(d==2){
      df=r*(p[1]+p[2]-r)+x2-pzeros
      #   BIC=log(dev_final)+log(n)*(r*(p[1]+p[2]-r)+x2)
  }else{
      df=r*(sum(p[1:(d-1)])-d+1)+x2-pzeros
      #  BIC=log(dev_final)+log(n)*(r*(sum(p[1:(d-1)])-d+1)+x2);
  }
  BIC=n*(log(2*pi)+log(dev_final))+n+log(n)*df
  
  tunning_index=list(deviance=n*(log(2*pi)+log(dev_final))+n,df=df)
  
  if(family=='binomial'){
    pre_error=sum(abs(as.numeric(linearpart_final>0)-y))/n
  }else{
    #pre_error is distance here
    pre_error=norm(as.matrix(y)-linearpart_final,'f')^2/sum(y^2)
  }
  rm(Md)
  #rm(Xj_eta)
  if(r>1){
    #rm(Xj_d)
  }
  
  if(Convergence_Index!=1){
    print("do not reach the stopping rule")
  }
  
  initialquantity_for_f=var(ob_f_differentinitial/n)
  initialquantity_for_loss_penalty=var(ob_f_loss_penalty_differentinitial/n)
  
  early_stop_for_lambda=0
  
  test=list(beta0=beta0_final,beta=beta_final,dev=dev_final,BIC=BIC,error=pre_error,
            num=num_final,ob_f_final=ob_f_final,ob_f_loss_penalty_final=ob_f_loss_penalty_final,
            ob_f_differentinitial=ob_f_differentinitial,
            ob_f_loss_penalty_differentinitial=ob_f_loss_penalty_differentinitial,
            linearpart_differentinitial=linearpart_differentinitial,
            initialquantity_for_f=initialquantity_for_f,
            initialquantity_for_loss_penalty=initialquantity_for_loss_penalty,tunning_index=tunning_index,early_stop_for_lambda=early_stop_for_lambda)
  return(test)
}













