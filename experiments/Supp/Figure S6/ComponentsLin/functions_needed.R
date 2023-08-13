#there are difference 
#for Tensor on  tensod
library(rTensor)
library(fda)
library(MASS)
library(glmnet)
library(cvTools)
library(splines2)
library(mpath)
library(restriktor)
library(quadprog)
#library(optiSolve)


Xj_function=function(X_sample,jj){
  Memorysave=0
  p1=dim(X_sample)
  p=p1[-length(p1)]
  d=length(p)
  if(class(X_sample)[1]=="Tensor"){
    TM=X_sample
  }else{
    TM=as.tensor(X_sample)
  } 
  
  for (dd in 1:d){
    Md[[dd]]=unfold(TM,c(d+1,dd),c(1:d)[-dd])
  }
  
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
        Xj=matrix(Md[[j]]@data%*%khatri_rao(beta[[d]],cumkr),n,p[j]*r)
      }
    }else{
      if(Memorysave==1){
        Md[[j]]=unfold(TM,c(d+1,j),c(1:d)[-j]) 
        Xj=matrix(Md[[j]]@data%*%khatri_rao(khatri_rao_list(beta[c(d:(j+1))]),cumkr),n,p[j]*r)
        Md[[j]]=0
        gc()
      }else{
        Xj=matrix(Md[[j]]@data%*%khatri_rao(khatri_rao_list(beta[c(d:(j+1))]),cumkr),n,p[j]*r)}
    }
  }
  
  if(jj==j)
  {
    return(Xj)
  }
  
  cumkr=khatri_rao(beta[[j]],cumkr)
  }
}



linearpart_penalty_function=function(y,X,beta,beta0,lambda){
X=X@data
D=full_R(beta)
p1=dim(X)
n=p1[length(p1)]
linearpart=c()
for(i in 1:n){
  linearpart[i]=sum(c(D)*c(as.array(X[,,,i])))+beta0[1]
}

l=loglike(y,linearpart,family = "gaussian")
P=l+n*lambda*penaltyvalue(beta,alpha=0,manifold=0)
L=l
test=list(L=L,lossP=P)
return(test)
}


loglike<-function(y,linearpart,family="binomial"){
  if(family=="binomial"){
    n=length(linearpart)
    dev_sum=0
    for(i in 1:n){
      if(linearpart[i]<=600){
        dev_sum=dev_sum+(log(1+exp(linearpart[i])))}else{
          dev_sum=dev_sum+linearpart[i]
        }
      
    }
    
    dev=t(as.matrix(y))%*%(as.matrix(linearpart))-sum(linearpart)-sum(dev_sum)
    return(dev)
  }
  #for guassian, it is equivalent to L2 loss
  if(family=="gaussian"){
    return(norm(as.matrix(y)-linearpart,'f')^2)
  }
}


full_R<-function(Beta){
  d=length(Beta)
  p=matrix(0,d,1)
  r=dim(Beta[[1]])[2]
  for(i in 1:d){
    p[i]=dim(Beta[[i]])[1]
  }
  B=array(0,c(p))
  for(j in 1:r){
    B1=Beta[[1]][,j]
    for(i in 2:d){
      B1=B1%o%Beta[[i]][,j]
    }
    B=B+B1
  }
  return(B)
}

arrange<-function(beta,a,manifoldfix=1){
  #the old version r and d are inter-changed.
  #r here is D+1
  ##r=length(beta)
  if(manifoldfix==1){
    d=length(beta)-1 
  }else{
  d=length(beta)
  }
  # d here is r
  ##d=dim(beta[[1]])[2]
  r=dim(beta[[1]])[2]
  if(a>=0){
    a=a^(1/d)
    for( i in 1:d){
      beta[[i]]=beta[[i]]*a
    }
  }else{
    a=(-a)^(1/d)
    for(i in 1:d){
      if(i==1){
        beta[[i]]=-beta[[i]]*a
      }else{
        beta[[i]]=beta[[i]]*a 
      }
    }
    
  }
  
  return(beta)  
}

#the rescale factor for elastic
#####this is used for Newton method in rescale strategy####
#f_d^\prime(x)
elastic_gradient_d=function(betad,alpha,lambda){
  betadl1=sum(abs(betad))
  betadl2square=sum((betad)^2) 
  
  numerator=0.5*4*(1-alpha)*betadl2square
  denominator1=-alpha*betadl1+sqrt(alpha^2*betadl1^2+4*lambda*(1-alpha)*betadl2square)
#  print(denominator1)
#  print("this")
#  print(betadl1)
#  print(lambda)
#  print(betadl2square)
#  print("that")
  denominator2=sqrt(alpha^2*betadl1^2+4*lambda*(1-alpha)*betadl2square)
  denominator=denominator1*denominator2
  
  gradient_for_d=numerator/denominator
  return(gradient_for_d)
  
}

#f_d(x)
elastic_f_lambda_d=function(betad,alpha,lambda){
  betadl1=sum(abs(betad))
  betadl2square=sum((betad)^2) 
  term1=log(-alpha*betadl1+sqrt(alpha^2*betadl1^2+4*lambda*(1-alpha)*betadl2square))
  term2=log(betadl2square)
  termd=term1-term2
  return(termd)
}

#f/f^\prime
elastic_foverfprime=function(beta,alpha,lambda,d,rr){
  numerator_f=0
  denominator_f_prime=0
  for(i in 1:d){
    betad=beta[[i]][,rr]
    numerator_f=numerator_f+elastic_f_lambda_d(betad,alpha,lambda)
    denominator_f_prime=denominator_f_prime+elastic_gradient_d(betad,alpha,lambda)
  } 
  numerator_f=numerator_f-d*log(2-2*alpha)
  foverfprime=numerator_f/denominator_f_prime
  return(c(foverfprime,numerator_f))
}

#Newton's method for r
#return the lagrange factor
elastic_newton_r=function(beta,alpha,d,rr){
 lambda=1
  for(i in 1:d){
    betad=beta[[i]][,rr]
   lambda=lambda*sum(betad^2)
  }
 lambda=lambda^(1/d)
 #the total iteration for newton is 10
 for(iter in 1:10){
   lambda=lambda-elastic_foverfprime(beta,alpha,lambda,d,rr)[1]
   #print("hehe")
   #print(elastic_foverfprime(beta,alpha,lambda,d,rr)[2])
   #later on add a if condition
   #
 }
 return(lambda)
}


#here lambda is the lagrange, not the tunning parameter in the model
#return the rescale factor by Newton, which need to be adjusted.
elastic_rescale_factor_rough=function(betad,alpha,lambda){
  betadl1=sum(abs(betad))
  betadl2square=sum((betad)^2)
  #the rescale factor
  numerator=-alpha*betadl1+sqrt(alpha^2*betadl1^2+4*lambda*(1-alpha)*betadl2square)
  denominator=2*(1-alpha)*betadl2square
  epsilond=numerator/denominator
  return(epsilond)
}


#return  the rescale factor
elastic_rescale_factor=function(beta,alpha,d,rr){
  epsilon_all=1
  lambda=elastic_newton_r(beta,alpha,d,rr)
 # print('test')
#  print(elastic_foverfprime(beta,alpha,lambda,d,rr)[2])
#  print('test')
  
  epsilon_d=c()
  for(i in 1:d){
    betad=beta[[i]][,rr]
    epsilon_d[i]=elastic_rescale_factor_rough(betad,alpha,lambda)
    epsilon_all=epsilon_all*epsilon_d[i]
  }
  #print(epsilon_all)
  #print(elastic_foverfprime(beta,alpha,lambda,d,rr)[2])
  tuningto1=(1/epsilon_all)^(1/d)
 #print(tuningto1)
 # print(epsilon_d)
  for(i in 1:d){
    epsilon_d[i]=epsilon_d[i]*tuningto1
  }
 # print(epsilon_d)
  return(epsilon_d)
}



#if beyond the limitation, then, no scale
#all=1 means rescale all components, otherwise, rescale the former d-1 components.
#l2=1 means
rescale<-function(beta,rescale_all=1,rescale_l2=1,elastic=0,alpha=0){
  if(rescale_all==1){
  d=length(beta)}else{
    #only rescale CP coefficients
    d=length(beta)-1
  }
  r=dim(beta[[1]])[2]
  lambda_r=rep(1,r)
  #in some step, part of beta has been chanted 
  #elastic rescale strategy is also included in the model
  beta_org=beta
  
  for(i in 1:r){
    for(ddd in 1:d ){
      if(rescale_l2==1){
      lambda_ir=sum(beta[[ddd]][,i]^2)^(0.5)}else{
      lambda_ir=sum(abs(beta[[ddd]][,i]))
      }
      #to avoid numerical problem
      if(is.nan(lambda_ir)==1){
        return(beta) 
      }
      #to avoid numerical problem
      if(lambda_ir==0){
        return(beta_org)
      }
      lambda_r[i]=lambda_r[i]*lambda_ir
      beta[[ddd]][,i]=beta[[ddd]][,i]/lambda_ir
    }
  }
  lambda_r=lambda_r^{1/d}
  for(i in 1:r){
    for(ddd in 1:d){
      beta[[ddd]][,i]=beta[[ddd]][,i]*lambda_r[i] 
    }
  }
  
  
  beta_elsatic=beta
  if(elastic!=0&&alpha!=1){
    for(rr in 1:r){
      rescaleforr=elastic_rescale_factor(beta,alpha,d,rr)
      for(ddd in 1:d){
        beta_elsatic[[ddd]][,rr]= beta[[ddd]][,rr]*rescaleforr[ddd]
        }
    
    }
    
    #to avoid numeraical problems
    if(is.na(sum(rescaleforr))!=1){
      #print("OK")
      #print(rescaleforr)
      beta=beta_elsatic
    }
  }
  

    
    
  
  return(beta)
}




#for the special case, i.e, X is a matrix,
#it will be changed
#note that we have changed the knots
BX_mat<-function(X,order=4,knots,type="truncated"){
  #note that we adjust the usual truncated power basis
  
  if(type=="bspline"){
    K=length(knots)+order-2
    p=dim(X)
    d=length(p)
    BX_array<-array(0,c(p,K))
    if(d==3){
      for(i_1 in 1:p[1]){
        for(i_2 in 1:p[2]){
          for(i_3 in 1:p[3]){
            BX_array[i_1,i_2,i_3,]<-bsplineS(X[i_1,i_2,i_3],knots,order,0)
          }
        }
      }}
    if(d==2){
      for(i_1 in 1:p[1]){
        for(i_2 in 1:p[2]){
          
          BX_array[i_1,i_2,]<-bsplineS(X[i_1,i_2],knots,order,0)
          
        }
      }
      
    }
    
    return(BX_array)
  }
  if(type=="truncated"){
    #K must be bigger than 4
    K=length(knots)+3-2
    P=dim(X)
    
    BX_array<-array(0,c(P,K))
    d=length(P)
    if(d==3){
      p1<-dim(X)[1]
      p2<-dim(X)[2]
      p3<-dim(X)[3]
      
      for(i_1 in 1:p1){
        for(i_2 in 1:p2){
          for(i_3 in 1:p3){
            Tr_basis=c()
            Tr_basis[1:3]=c(X[i_1,i_2,i_3],X[i_1,i_2,i_3]^2,X[i_1,i_2,i_3]^3)
            for(i in 4:K){
              if(X[i_1,i_2,i_3]-knots[i-3]>0){
                Tr_basis[i]=(X[i_1,i_2,i_3]-knots[i-3])^3
              }
              if(X[i_1,i_2,i_3]-knots[i-3]<=0){
                Tr_basis[i]=0
              }
            }
            BX_array[i_1,i_2,i_3,]<-Tr_basis
          }
          
        }
        
      }
    }
    if(d==2){
      p1<-dim(X)[1]
      p2<-dim(X)[2]
      for(i_1 in 1:p1){
        for(i_2 in 1:p2){
          
          Tr_basis=c()
          Tr_basis[1:3]=c(X[i_1,i_2],X[i_1,i_2]^2,X[i_1,i_2]^3)
          for(i in 4:K){
            if(X[i_1,i_2]-knots[i-3]>0){
              Tr_basis[i]=(X[i_1,i_2]-knots[i-3])^3
            }
            if(X[i_1,i_2]-knots[i-3]<=0){
              Tr_basis[i]=0
            }
          }
          BX_array[i_1,i_2,]<-Tr_basis
          
          
        }
        
      }
      
    }
    
    return(BX_array)  
    
    
  }
}

#Suppose the sample is a 3-order array
#actually, the sample array works too, but I want to make a consistent form, then, then tensor format is used
BX_sample_fun<-function(X_sample,K=6,order=4,type="truncated"){
  
  
  X_vector=as.vector(X_sample)
  knots=quantile(X_vector,probs=c(seq(0,1,1/(K-3+1))))
  p=dim(X_sample)
  d=length(p)-1
  n=p[d+1]
  if(type=="bspline"){
    K=length(knots)+order-2
  }
  BX_sample_array<-array(0,c(p[-length(p)],K,n))
  if(d==3){
    for(i in 1:n){
      BX_sample_array[,,,,i]<-BX_mat(X_sample[,,,i],order,knots,type=type)
    }}
  if(d==2){
    for(i in 1:n){
      BX_sample_array[,,,i]<-BX_mat(X_sample[,,i],order,knots,type=type)
    }} 
  
  
  BX_sample_tensor<-as.tensor(BX_sample_array)
  knots_and_BX_sample_tensor<-list(knots=knots,BX_sample_tensor=BX_sample_tensor)
  return(knots_and_BX_sample_tensor)
}

prediction_function_nonlinear<-function(b0,BB,X,y,family='gaussian',hat=0){
  d=length(dim(X))-1
  if( family=='gaussian'){
    n=length(y)
    BB=as.tensor(BB)
    if(class(X)[1]=="Tensor"){
    }else{
      X=as.tensor(X)
    }
    yhat=c()
    for(i in 1:n){
      if(d==4){
      yhat[i]=innerProd(BB,X[,,,,i])+b0
      }
      if(d==3){
        yhat[i]=innerProd(BB,X[,,,i])+b0  
      }
    }
    rm(X)
    gc()
    if(hat==0){
    MSE=norm(as.matrix(y)-yhat,'f')^2/n
    return(MSE)
    }else{
        MSE=norm(as.matrix(y)-yhat,'f')^2/n
        test=list(MSE=MSE,yhat=yhat)
        return(test)
    }
    
  }
}
prediction_function_linear<-function(b0,BB,X,y,family='gaussian'){
  if( family=='gaussian'){
    n=length(y)
    if(class(BB)[1]=="Tensor"){
    }else{
      BB=as.tensor(BB)
    }
    
    
    if(class(X)[1]=="Tensor"){
    }else{
      X=as.tensor(X)
    }
    
    yhat=c()
    d=length(dim(X))-1
    if(d==3){
    for(i in 1:n){
      yhat[i]=innerProd(BB,X[,,,i])+b0
    }
    }
    if(d==2){
    for(i in 1:n){
        yhat[i]=innerProd(BB,X[,,i])+b0
        }
    }
    #RPE
    MSE=norm(as.matrix(y)-yhat,'f')^2/n
    return(MSE)
  }
}

rlm<-function(X,y,uk)
{
  return(0)
}
#X_sample=array(rnorm(10*5*6*20),c(10,5,6,20))


#optimization on the oblique manifold

grad<-function(Y,Xk_vec,D,K,r){
  
  nablaLvec=(-2*t(D)%*%Y+2*t(D)%*%D%*%Xk_vec)
  nablaL=matrix(nablaLvec,K,r)
  Xk=matrix(Xk_vec,K,r)
  
  #print(dim(Xk))
  #print(dim(nablaL))
  if(r==1){
    gradf=nablaL-Xk%*%(t(Xk)%*%nablaL) 
  }else{
  gradf=nablaL-Xk%*% diag(diag(t(Xk)%*%nablaL))
  }
  return(gradf)
}
Retraction<-function(grad,alpha,Xk){
Xk1=Xk+alpha*(-grad)
r=dim(Xk1)[2]
if(r==1){
  Xk1=Xk1%*%(t(Xk1)%*%Xk1)^(-0.5) 
}else{
Xk1=Xk1%*%diag((diag(t(Xk1)%*%Xk1))^(-0.5))
}
return(Xk1)
}

falphaL=function(D,X,Y){
  return(sum((Y-D%*%c(X))^2))
}

update_alpha=function(fast=0,alphafast=0,Y,BD1,D,K,r,maxupdate_out=1,maxupdate_in=40,scalefactor=0.5,thresholdformani=1e-10,guess=0){
  X0=BD1
  #print(dim(D))
  #print(dim(X0))
  #print(dim(Y))
  f0=falphaL(D,X0,Y)
  Xl=X0
  fl=f0
  fcontrol=f0

  gradX0=grad(Y,c(X0),D,K,r)
  gradXl=gradX0
  
  #the outer loop, both for fast or not fast
  for(l  in 1:maxupdate_out){
  #print(l)  
  norm2gradXl=norm(gradXl,'f')^2
  #print(norm2gradXl)
  
  #the basic method to choose alpha, guess or 0.5^l*alpha
  if(guess==1){
    alpha=1/(norm2gradXl^0.5)
    print(alpha)
    print(fl)
  }else{
  alpha=2
  }
  
  
  #this will be used for those cases that do not satisfy Armijo
  Xmay=Xl
  fmay=fl
  #this is used for those cases that the monifold dose not move
  #fcontrol=fl
  
  #index for Armijo criterion
  ArmijoControl=0
  
  #The fast setting
  ArmijoControlfast=0
  alphafastuse=alphafast
  
  #this two indexes are used to 
  controlindexmore=0
  controlindexless=0
  
  #The fast setting
  if(fast==1){
  #this loop is to find a small step size
  for(i in 1:maxupdate_in){
  Xl1=Retraction(grad=gradXl,alpha=alphafastuse,Xk=Xl) 
  fl1=falphaL(D=D,X=Xl1,Y=Y)
  
  if((fl1-fl)<=-0.5*alphafastuse*norm2gradXl){
    if(controlindexless==1){
      #print("Armijofast")
      ArmijoControlfast=1
      Xl=Xl1
      fl=fl1
      gradXl=grad(Y,c(Xl),D,K,r)
      break;
    }
    controlindexmore=1
    alphafastuse=2*alphafastuse
    
  }else{
    controlindexless=1
    alphafastuse=0.5*alphafastuse
    
  }
  }
    }
  
  if(ArmijoControlfast==1){
    alphafast=alphafastuse
    break;
  }
  
  
  #the old and basic one
  for(i in 1:maxupdate_in){
    if(guess==1){
      #this guess may need to be adjusted
      alpha=alpha*scalefactor
    }else{
      alpha=alpha*scalefactor
    } 
     

  #print(alpha)
  Xl1=Retraction(grad=gradXl,alpha=alpha,Xk=Xl)
  #print(Xl1)
  #print(dim(D))
  #print(dim(Xl1))
  fl1=falphaL(D=D,X=Xl1,Y=Y)
  #print(fl1)
  if((fl1-fl)<=(-0.5*alpha*norm2gradXl)){
    ArmijoControl=1
    Xl=Xl1
    fl=fl1
    #print(fl1)
    alphafast=alpha
    gradXl=grad(Y,c(Xl),D,K,r)
    break;}
  
  #if decline, we can accept it
   if(fl1<fl){
     Xmay=Xl1
     fmay=fl1
     alphafast=alpha
   }
  
  }
  
  
  #print(fl1-fl)
  #print(-0.5*alpha*norm2gradXl)
  
  
  #not use Armijo step-size rule 
  if(ArmijoControl==1){
#print("Armijo")
  }else{
  if(fmay>=fcontrol){

#    print("manifold optimization cannot move")
    Notmove=list(Xl=X0,alphafast=1)
    return(Notmove)
    break;
    }else{
#print("Not Armijo")
#print(fl)
#print(fmay)
    Xl=Xmay
    fl=fmay
    gradXl=grad(Y,c(Xl),D,K,r)
    }
#    
  }
#  
 
  if(max(abs(gradXl))<=thresholdformani*(1+max(abs(gradX0)))){
    print('good')
     break;}
 
  }
 
  if(fl>f0){
    result=list(Xl=X0,alphafast=alphafast,fl=f0) 
  }else{
  result=list(Xl=Xl,alphafast=alphafast,fl=fl)
  }
  return(result)
}


#only used when l2 penalty
#later, may works for other penalty. So we can just keep alpha at present.
penaltyvalue=function(beta,alpha,manifold=1,QCQP=0){
  if(manifold==1||QCQP==1){
  D=length(beta)-1
  }else{
    D=length(beta)
  }
  
 
  value=0
  value1=0
  value2=0
  if(alpha==0){
  for(i in 1:D){
    value=value+sum((beta[[i]])^2)
  }
  }else{
    for(i in 1:D){
      value1=value1+sum(abs(beta[[i]]))
    }
    value1=2*value1
    for(i in 1:D){
        value2=value2+sum((beta[[i]])^2)
    }
    value2=value2
    
    value=value1+value2
    
  }
  return(value)
}


distancefunction=function(linearpartlist){
  d=length(linearpartlist)
  distancematrix=matrix(0,d,d)
  for(i in 1:d){
    for(j in i:d){
      distancematrix[i,j]=(sum((linearpartlist[[i]]-linearpartlist[[j]])^2))^0.5
    }
  }
  return(distancematrix)
}

downsizearray=function(Xarray,pobject){
  p=dim(Xarray)
  D=length(p)
  d=length(pobject)
  Xarray_new=array(0,pobject)
  index=list()
  for(i in  1:d){
    common_difference=(p[i]-1)/(pobject[i]-1)
    index[[i]]=floor(seq(1,p[i],common_difference))
  }
  
  #assume D=d now. in the future, this may be improved.
  if(d==2){
    Xarray_new=Xarray[c(index[[1]]),c(index[[2]])]
  }
  if(d==3){
    Xarray_new=Xarray[c(index[[1]]),c(index[[2]]),c(index[[3]])]
  }
  if(d==4){
    Xarray_new=Xarray[c(index[[1]]),c(index[[2]]),c(index[[3]]),c(index[[4]])] 
  }
  test=list(Xarray_new=Xarray_new,index=index)
  return(test)
}

uparraysize=function(Xarray_new,pobjectup,index){
  Xarray=array(0,pobjectup)
  d=length(pobjectup)
  
  if(d==2){
    Xarray[c(index[[1]]),c(index[[2]])]=Xarray_new
  }
  if(d==3){
    Xarray[c(index[[1]]),c(index[[2]]),c(index[[3]])]=Xarray_new
  }
  if(d==4){
    Xarray[c(index[[1]]),c(index[[2]]),c(index[[3]]),c(index[[4]])] =Xarray_new
  }
  return(Xarray)
}

#This is following hua zhou's code
linspace=function(d1,d2,n){
 n=floor(n) 
 n1=n-1
 c=(d2-d1)*(n1-1)
 if(c==Inf){
   if((d2-d1)==Inf){
     y=d1+(d2/n1)*c(0:n1)-d1/n1*(0:n1)
   }else{
     y=d1+c(0:n1)*((d2-d1)/n1)
   }
 }else{
     y=d1+c(0:n1)*(d2-d1)/n1
 }
 if(is.na(y[1])){
 }else{
   if(d1==d2){
     y=d1
   }else{
     y[1]=d1
     y[length(y)]=d2
   }
 }
 return(y)
 }
 
resize_array=function(Xarray,pobject){
  p=dim(Xarray)
  D=length(p)
  d=length(pobject)
  U=list()
  Xarray_new=Xarray
  if(class(Xarray_new)[1]=="Tensor"){
    tensor_indicator=1
  }else{
    tensor_indicator=0
  Xarray_new=as.tensor(Xarray_new)  
  }
  for(d in 1:D){
    xi=linspace(1,p[d],pobject[d])
    j1=floor(xi)
    j2=j1+1
    w1=j2-xi
    w2=1-w1
    j2[length(j2)]=p[d]
    Ud=as.matrix(spMatrix(pobject[d],p[d],c(1:pobject[d],1:pobject[d]),c(j1,j2),c(w1,w2)))
    Xarray_new=ttm(Xarray_new,Ud,m =d)
  }
  #return the same class as Xarray
  if(tensor_indicator==1){
    return(Xarray_new)
  }else{
      return(Xarray_new@data)
    }
}


#someting new 
dowsizearray_blocker=function(X_array,sizeofslinding=NA){
  if(class(X_array)[1]=="Tensor"){
    tensor_indi=1
  }else{
    tensor_indi=0
  }
  
  p1=dim(X_array)
  D1=length(p1)
  n=p1[D1]
  K=p1[D1-1]
  p=p1[1:(D1-2)]
  p_aug=p
  D=length(p)
  #Up to now, only support D=2, D=3. It is easy to extended to D=4,5...
  if(is.na(sizeofslinding)){
  if(D==2){
  if(p[1]%%2!=0){
  p_aug[1]=p[1]+1
  }
  if(p[2]%%2!=0){
  p_aug[2]=p[2]+1  
  }
    
  X_array_middle=array(0,c(p_aug,K,n))
  #print(dim(X_array))
  #print(dim(X_array_middle))
  if(tensor_indi==1){
    X_array_middle[1:p[1],1:p[2],1:K,1:n]=X_array@data 
  }else{
  X_array_middle[1:p[1],1:p[2],1:K,1:n]=X_array
  }
  X_array= X_array_middle
  X_array_new=array(0,c(0.5*p_aug[1],0.5*p_aug[2],K,n))  
  for(i in 1:(0.5*p_aug[1])){
    for(j in 1:(0.5*p_aug[2])){
      X_array_new[i,j,,]=0.25*(X_array[(2*(i-1)+1),(2*(j-1)+1),,]+X_array[(2*(i-1)+1),(2*j),,]+
                                 X_array[(2*i),(2*(j-1)+1),,]+X_array[(2*i),(2*j),,])
    }
  }  
  }
  
  if(D==3){
   if(p[1]%%2!=0){
      p_aug[1]=p[1]+1
    }
    if(p[2]%%2!=0){
      p_aug[2]=p[2]+1  
    }
    if(p[3]%%2!=0){
      p_aug[3]=p[3]+1  
    }
    
    X_array_middle=array(0,c(p_aug,K,n)) 
    if(tensor_indi==1){
      X_array_middle[1:p[1],1:p[2],1:p[3],1:K,1:n]=X_array@data 
    }else{
      X_array_middle[1:p[1],1:p[2],1:p[3],1:K,1:n]=X_array
    }
    X_array= X_array_middle
    X_array_new=array(0,c(0.5*p_aug[1],0.5*p_aug[2],K,n))  
    for(i in 1:(0.5*p_aug[1])){
      for(j in 1:(0.5*p_aug[2])){
        for(k in 1:(0.5*p_aug[3])){
        X_array_new[i,j,k,,]=0.125*(X_array[(2*(i-1)+1),(2*(j-1)+1),(2*(k-1)+1),,]+X_array[(2*(i-1)+1),(2*j),(2*(k-1)+1),,]+
                                      X_array[(2*i),(2*(j-1)+1),(2*(k-1)+1),,]+X_array[(2*i),(2*j),(2*(k-1)+1),,]+
                                      X_array[(2*(i-1)+1),(2*(j-1)+1),(2*(k-1)+2),,]+X_array[(2*(i-1)+1),(2*j),(2*(k-1)+2),,]+
                                      X_array[(2*i),(2*(j-1)+1),(2*(k-1)+2),,]+X_array[(2*i),(2*j),(2*(k-1)+2),,])
        }
      }
    }    
    
    
    
  }    
    
  }
  if(tensor_indi==1){
    X_array_new=as.tensor(X_array_new)
  }
  return(X_array_new)
}

upsizearray_blocker=function(B_reduce,sizeofslinding=NA){
p=dim(B_reduce)
D=length(p)
B_aug=array(0,c(2*p[-D],p[D]))

if(is.na(sizeofslinding)){
if(D==3){
  for(i in 1:p[1]){
    for(j in 1:p[2]){
      B_aug[(2*(i-1)+1):(2*i),(2*(j-1)+1):(2*j),]=B_reduce[i,j,]
    }
  }
  
  

}

if(D==4){
for(i in 1:p[1]){
  for(j in 1:p[2]){
    for(k in 1:p[3]){
      B_aug[(2*(i-1)+1):(2*i),(2*(j-1)+1):(2*j),(2*(k-1)+1):(2*k),]=B_reduce[i,j,k,] 
    }
  }
} 
  
}    
  
}
  
return(B_aug)  
}


coodinateforalpha<-function(y,DX,alphaold,alphaTolFun=0.01,Maxiterforalpha=100000){
  r=dim(alphaold)[2]
  K=dim(alphaold)[1]
  alpha=alphaold
  dev0=sum((y-DX%*%c(alphaold))^2)
  
  
  for (iter in 1:Maxiterforalpha) {

    
    for(i in 1:r){
      #Qr_vec=rep(0,r*K)
      #Qr_vec[((i-1)*K+1):(i*K)]=1
      #Qr=diag(Qr_vec)
      
      Qr=diag(K)
      C_Qr=quadcon(Qr,val=1)
      y_new=y-DX[,-(((i-1)*K+1):(i*K))]%*%c(alpha[,-i]) 
      DXX=DX[,((i-1)*K+1):(i*K)]
   
      
      goal=sum((y-DX%*%c(alpha))^2)
      print(goal) 
      
      goal2=sum((y_new-DXX%*%c(alpha[,i]))^2)
      print(goal2)      
      for(j in 1:1){
        print("old")
        print(sum((y_new-DXX%*%c(alpha[,i]))^2))
        print(alpha[,i])
      fit=glmnet(DXX,y_new,lambda = 0.5,alpha=0,intercept = 0)
      #y_new=y_new-coef(fit)[1]
      alpha[,i]=coef(fit)[-1]
      print("new")
      print(sum(y_new^2))
      print(sum((y_new-DXX%*%c(alpha[,i]))^2))
      #print(alpha[,i])
      if(sum(alpha[,i]^2)<=10){
        print("haha")
        break;
      }
      }
      
      
      #goal=sum((y-DX%*%c(alpha))^2)
      #print(goal) 
      
      #DXXX=t(DXX)%*%DXX
      #Q_F=quadfun(Q=DXXX,a=-2*t(DXX)%*%y_new,d=sum(y_new^2))
      #print(c(alpha[,i])%*%DXXX%*%c(alpha[,i])-2*c(alpha[,i])%*%t(DXX)%*%y_new+sum(y_new^2))
      #mycop=cop(f=Q_F,qc=C_Qr)  
      #res <- solvecop(mycop,solver = "csdp",quiet=TRUE)
     # alpha[,i]=res$x
    } 
    
   
    devtmp=sum((y-DX%*%c(alpha))^2)
    devdiff=dev0-devtmp
    dev0=devtmp
    if(abs(devdiff)<alphaTolFun*(devtmp+1)){
      return(alpha)
    }  
  }

  
}


dualascentforalpha<-function(y,DX,alphaold,C=1,Maxiterfordual=1000,TolFundualascent=1e-7,MaxIterInner=40){
r=dim(alphaold)[2]
K=dim(alphaold)[1] 
initialvalue=c()
h_alpha_intial=c()
initialsacale=2
#initial value

for(i in 1:10){
  initialpointtest=ginv(t(DX)%*%DX+(initialsacale^i)*diag(r*K))%*%t(DX)%*%y 
  
  for(j in 1:r){
    h_alpha_intial[j]=sum(initialpointtest[((j-1)*K+1):(j*K)]^2)-C 
  }
    initialvalue[i]=abs(sum(initialpointtest*initialpointtest)-r*C)
   initialindex=sum(h_alpha_intial<=0)
   if(initialindex==r){
    
   }else{
     initialvalue[i]=Inf
   }
    
  if(initialvalue[i]<0.1){
    break
  }
}
initial_index=which.min(initialvalue)
#print(initial_index)

lambda0=rep(initialsacale^initial_index,r)
lambda=lambda0

n=dim(DX)[1]

lambda_diag_vec=rep(1,r*K)
for(i in 1:r){
lambda_diag_vec[((i-1)*K+1):(i*K)]=lambda[i]*lambda_diag_vec[((i-1)*K+1):(i*K)]
}
lambda_diag=diag(lambda_diag_vec)
h_alpha=rep(0,r)


test_lambda=list()
test_L_x=c()
test_f=c()
test_alpha=list()

test_halpha=list()

for(iter in 1:Maxiterfordual){

if(iter==1){  
alphavec=ginv(t(DX)%*%DX+lambda_diag)%*%t(DX)%*%y
for(i in 1:r){
h_alpha[i]=sum(alphavec[((i-1)*K+1):(i*K)]^2)-C 
}

f_0=sum((y-DX%*%alphavec)^2)
L_x_lambda_0=f_0+sum(h_alpha*lambda)

test_lambda[[1]]=lambda
test_alpha[[1]]=alphavec
test_L_x[1]=L_x_lambda_0

test_f[1]=f_0

test_halpha[[1]]=h_alpha

}else{
  tk=2
  L_x_lambda_tem_final=L_x_lambda_0
  indexarmijo=0
  h_alpha_final=h_alpha
  lambda_final=lambda
  f0_final=f_0
  alphavec_final=alphavec
  
  for(j in 1:MaxIterInner){
  tk=tk*0.5
  
  for(i in 1:r){
  #h_alpha[i]=sum(alphavec[((i-1)*K+1):(i*K)]^2)-C
   lambda[i]=max((lambda[i]+tk*h_alpha[i]),0)
  }
  
  lambda_diag_vec=rep(1,r*K)
  for(i in 1:r){
  lambda_diag_vec[((i-1)*K+1):(i*K)]=lambda[i]*lambda_diag_vec[((i-1)*K+1):(i*K)]
    }
  
  lambda_diag=diag(lambda_diag_vec)
  
  alphavec=ginv(t(DX)%*%DX+lambda_diag)%*%t(DX)%*%y
  
  f_0=sum((y-DX%*%alphavec)^2)
  
  for(i in 1:r){
    h_alpha[i]=sum(alphavec[((i-1)*K+1):(i*K)]^2)-C
    #lambda[i]=max((lambda[i]+tk*h_alpha[i]),0)
  }
  
  L_x_lambda_tem_before=f_0+sum(h_alpha*lambda)
  if(L_x_lambda_tem_before-L_x_lambda_0>(tk*0.5*sum(h_alpha^2))){
    #feasbile  adjustment
    if(sum(h_alpha<=0)==r){
    indexarmijo=1
    L_x_lambda_tem=L_x_lambda_tem_before
  
    break;
    }else{}
  }
  
  if(L_x_lambda_tem_before>L_x_lambda_tem_final){
    L_x_lambda_tem_final=L_x_lambda_tem_before
    h_alpha_final=h_alpha
    lambda_final=lambda
    f0_final=f_0
    alphavec_final=alphavec
  }
  
  }
  
  if(indexarmijo==1){
    #print("hehe")
  }else{
    h_alpha=h_alpha_final
    lambda=lambda_final
    L_x_lambda_tem=L_x_lambda_tem_final
    f_0=f0_final
    alphavec=alphavec_final
  }
    
  
  L_x_lambda_diff=L_x_lambda_tem-L_x_lambda_0
  L_x_lambda_0=L_x_lambda_tem

  
  
  
  test_lambda[[iter]]=lambda
  test_alpha[[iter]]=alphavec

  
  test_halpha[[iter]]=h_alpha
  
  feasibleindex=sum(h_alpha<=(0))
  if(feasibleindex==r){
  test_L_x[iter]=L_x_lambda_0
  test_f[iter]=f_0
  }else{
    test_L_x[iter]=-Inf 
    test_f[iter]=Inf
  }
  
  if(max(lambda)==0){
    break;
  }
  
  if(abs(L_x_lambda_diff)<=TolFundualascent*(L_x_lambda_0)){
    if(iter>=10){
    break;
      }else{}
  }

}
}



indexf=which.max(test_L_x)

#indexf=which.min(test_f)
index=indexf
if(sum(test_halpha[[indexf]]<=0)==r){

alphavec_final=test_alpha[[indexf]]
}else{
  
  alphavec_final=c(alphaold)
  #alphavec_final=test_alpha[[indexf]]
}

res=list(test_alpha=test_alpha,alphavec_final=alphavec_final,test_L_x=test_L_x,test_f=test_f,test_lambda=test_lambda,test_halpha=test_halpha,index=index)

return(res)
}


cv_kruskal_sparsereg=function(X_sample,y,fold=5,alpha,lambdasequence=c(0,0.0001),r=2,initialpoints=5,warmstart=0,manifold=0,QCQP=0,TolFun=0.001){
 n=length(y)
 num_of_tunning=length(lambdasequence)
 set.seed(2019)
 idgroup=split(sample(n,n,replace = FALSE),rep(1:fold, length = n))
 MSE_train=matrix(0,fold,num_of_tunning)
 MSE_test=matrix(0,fold,num_of_tunning)
 lambdasequence=c(0,0.00001,0.0001,0.0005,0.0008,0.001,0.002,0.01,0.05)
 
 if(alpha==0){
 rescalparameter=1
 }else{
 rescalparameter=0 
 }
 
 for(i in 1:fold){
   idtest=idgroup[[i]]
   y_train=y[-idtest]
   y_test=y[idtest]
   X_sample_train=X_sample[,,,-idtest]
   X_sample_test=X_sample[,,,idtest]
  
   for(iter in 1:num_of_tunning){
     
     fit_nonlinear_list=broadcasted_sparsetenreg(TolFun=TolFun,rescale_all=0,rescale_l2=rescalparameter,r=r,X_sample_train,y_train,penalty=c("L1"),knots=knots,restriction = 0,lambda=lambdasequence[iter],gamma=0,alpha_gamma=0,
                                                            alpha=alpha,Replicates = initialpoints,family="gaussian",MaxIter = 1000,manifold=manifold,warmstart=warmstart,QCQP = QCQP)
     
     BB=full_R(fit_nonlinear_list$beta)
     b0=fit_nonlinear_list$beta0
     MSE_train[[i,iter]]=prediction_function_nonlinear(b0,BB,X_sample_train,y_train)
     MSE_test[[i,iter]]=prediction_function_nonlinear(b0,BB,X_sample_test,y_test)
   }
   
   
   
   

 }
 
 MSE_train_mean=apply(MSE_train,2,sum)/fold
 MSE_test_mean=apply(MSE_test,2,sum)/fold
 
 test=list(MSE_train=MSE_train,MSE_test=MSE_test,MSE_train_mean=MSE_train_mean,MSE_test_mean=MSE_test_mean,lambdasequence=lambdasequence)
 return(test)
  
}


#sum((y-DX%*%c(alphaold))^2)
#sum((y-DX%*%hehe$alphavec_final)^2)

#hehe=dualascentforalpha(y,DX,alphaold,C=1,Maxiterfordual=1000,TolFundualascent=0.00001)
#sum((y-DX%*%c(alphaold))^2)
#sum((y-DX%*%hehe$alphavec_final)^2)

#alpha_final=matrix(hehe$alphavec_final,10,4)
#sum(alpha_final[,1]^2)



#dada=glmnet(DX,y,alpha=0,lambda=10,family = "gaussian",intercept = 0)
#alphaglmnet=coef(dada)[-1]
#sum((y-DX%*%alphaglmnet)^2)
#alpha_final=matrix(alphaglmnet,10,4)
#sum(alpha_final[,2]^2)


#DX=matrix(rnorm(40000),1000,40)
#alphaold_vec=c(rnorm(40))
#alphaold=matrix(alphaold_vec,10,4)
#y=DX%*%alphaold_vec+rnorm(1000)
#sum(y^2)

#sum(alphaold[,3]^2)

#sum((y-DX%*%c(alphaold))^2)

#lambda0=100*rep(1,r)
#lambda0[2]=5*lambda0[2]
#lambda=lambda0

#n=dim(DX)[1]

#lambda_diag_vec=rep(1,r*K)
#for(i in 1:r){
#  lambda_diag_vec[((i-1)*K+1):(i*K)]=lambda[i]*lambda_diag_vec[((i-1)*K+1):(i*K)]
#}
#lambda_diag=diag(lambda_diag_vec)
#h_alpha=rep(0,r)

#alphavec=ginv(t(DX)%*%DX+1*lambda_diag)%*%t(DX)%*%y
#alphatest=matrix(alphavec,10,4)
#sum(alphatest[,1]^2)
#sum(alphatest[,2]^2)
#sum(alphatest[,3]^2)
#sum(alphatest[,4]^2)

Socpforalpha=function(y,DX,alphaold,C=1,Maxiterfordual=30,TolFundualascent=0.001,MaxIterInner=40){
r=dim(alphaold)[2]
K=dim(alphaold)[1]
n=length(y)
Ar_list=list()
g_list=list()
d_list=list()
f_vec=c()
for(i in 1:r){
Ar_vec=rep(0,r*K+1)
Ar_vec[((i-1)*K+1):(i*K)]=1
Ar_list[[i]]=diag(Ar_vec)[1:(r*K),]
g_list[[i]]=rep(0,r*K)
d_list[[i]]=rep(0,r*K+1)
f_vec[i]=C
}
#Corresponding the quadratic 
D_aug=matrix(0,n+1,r*K+1)
D_aug[2:(n+1),1:(r*K)]=DX
D_aug[1,r*K+1]=-0.5
Ar_list[[r+1]]=D_aug
g_aug=rep(0,n+1)
g_aug[1]=0.5
g_list[[r+1]]=g_aug
d_aug=rep(0,r*K+1)
d_aug[r*K+1]=0.5
d_list[[r+1]]=d_aug
f_vec[r+1]=0.5

c_cccp_list=list()
for(i in 1:(r+1)){
F1=Ar_list[[i]]
g1=g_list[[i]] 
d1=d_list[[i]]
f1=f_vec[i] 
c_cccp_list[[i]]=socc(F = F1, g = g1, d = d1, f = f1)    
}

q_use=c(-2*t(y)%*%DX,1)

ans <- cccp(q = q_use, cList = c_cccp_list, optctrl = ctrl())

res=getx(ans)

return(res)  

}

Scsforalpha<-function(y,DX,alphaold,C=1){
  obj=c(-2*t(y)%*%DX,1)
  
  r=dim(alphaold)[2]
  K=dim(alphaold)[1]
  n=length(y)
  Ar_list=list()
  g_list=list()
  d_list=list()
  f_vec=c()
  for(i in 1:r){
    Ar_vec=rep(0,r*K+1)
    Ar_vec[((i-1)*K+1):(i*K)]=1
    Ar_list[[i]]=diag(Ar_vec)[1:(r*K),]
    g_list[[i]]=rep(0,r*K)
    d_list[[i]]=rep(0,r*K+1)
    f_vec[i]=C
  }
  
  #Corresponding the quadratic 
  D_aug=matrix(0,n+1,r*K+1)
  D_aug[2:(n+1),1:(r*K)]=DX
  D_aug[1,r*K+1]=-0.5
  Ar_list[[r+1]]=D_aug
  g_aug=rep(0,n+1)
  g_aug[1]=0.5
  g_list[[r+1]]=g_aug
  d_aug=rep(0,r*K+1)
  d_aug[r*K+1]=0.5
  d_list[[r+1]]=d_aug
  f_vec[r+1]=0.5
  
  AA=matrix(0,((r*K+1)*r+n+2),r*K+1)
  bb=c()
  for(i in 1:r){
  AA[((i-1)*(K*r+1)+1):((i-1)*(K*r+1)+K*r),]=-Ar_list[[i]]
  AA[((i-1)*(K*r+1)+K*r+1),]=0
  
  bb[((i-1)*(K*r+1)+1):((i-1)*(K*r+1)+K*r)]=0
  bb[((i-1)*(K*r+1)+K*r+1)]=C
  }
  
  AA[(r*(K*r+1)+1):(r*(K*r+1)+n+1),]=-D_aug
  AA[(r*(K*r+1)+n+2),]=-d_aug
  
  bb[(r*(K*r+1)+1):(r*(K*r+1)+n+1)]=g_aug
  bb[(r*(K*r+1)+n+2)]=0.5

  b=bb
  A=AA
  cone<-list(q=c((r*K+1),(r*K+1),(n+2)))
  #cone=list(q=r*K+1+r*K+1+n+2)

  


  
  control <- list(eps = 1e-7, max_iters = 5000)
  sol<-scs(A,b,obj,cone)
  res=sol$x
  
  fl=sum((y-DX%*%res[-length(res)])^2)
  test=list(x=res,fl=fl)
  return(sol)
}

#library(cccp)
#library(scs)
#DX=matrix(rnorm(100),25,4)
#y=DX%*%c(3.3,1.5,-4.1,1.6)+rnorm(25)
#alphaold=matrix(rnorm(4),2,2)
#hehe=Socpforalpha(y,DX,alphaold,C=1)
#sum((y-DX%*%hehe[-5])^2)
#Scsforalpha(y,DX,alphaold,C=1)

#ee=dualascentforalpha(y,DX,alphaold)
#sum((y-DX%*%c(ee$alphavec_final))^2)

warmstartdowsize=function(warmstart=3,r=2,X_sample=NA,shrink_factor_scale=5,B0=NA,B0_list=NA,BB0=NA,BB0_list=NA){
#From a lot of  experiment, warmstart=3 is the best
  p1=dim(X_sample)
  n=p1[length(p1)]
  p=p1[-length(p1)]
  K=p[length(p)]
  d=length(p)
  
  
  B0_list_reduce=B0_list
  B0_reduce=B0
  BB0_reduce=BB0
  BB0_list_reduce=BB0_list
  
  X_sample_reduce=X_sample

  shrink_factor=(n/shrink_factor_scale)/(r*sum(p))

  if(warmstart==3){
      print(shrink_factor)
    if(shrink_factor<1){
#      print(shrink_factor)
#      print(p)
#      print(c(p[1:(d-1)]*shrink_factor,K))
#      heheda=c(p[1:(d-1)]*shrink_factor,K)
#      print(floor(heheda))
      #there is a bug  for R package itselt. To be more specific, it is  the floor function.
      p_reduce=floor(c(p[1:(d-1)]*shrink_factor,K))
      #if(p_reduce[1]<=1||p_reduce[2]<=1||p_reduce[3]<=1){
      #  p_reduce=floor(c(p[1:(d-1)]*shrink_factor+1,K))
      #}
      #      if(d==4){
      #   if(p_reduce[3]<=2){
      #        p_reduce[3]=p[3]
      #    }
      #    if(p_reduce[2]<=2){
      #       p_reduce[2]=3
      #    }
      #    if(p_reduce[1]<=2){
      #        p_reduce[1]=3
      #    }
      #}
      #if(d==5){
      #   for(test in 1:4)
      #    if(p_reduce[test]<=2){
      #        p_reduce[test]=3
      #    }
      #}
      for(test in 1:(d-1)){
          if(p_reduce[test]<=5){
            p_reduce[test]=min(p[test],3)
          }
      }
      
      
      
      

      
      #print(p_reduce)
      X_sample_reduce=resize_array(X_sample,c(p_reduce,n))
      r_reduce=r
      #initial parameter downsize  
      if(is.na(sum(B0_list[[1]][[1]]))!=1){
        Number_ini=length(B0_list)
        B0_list_reduce=B0_list
        for(rep in 1:Number_ini){
          for(sizepointer in 1:d){
            middle=B0_list[[rep]][[sizepointer]]
            print(dim(middle))
            print(p_reduce[sizepointer])
            B0_list_reduce[[rep]][[sizepointer]]=middle[1:p_reduce[sizepointer],1:r]
          }
        }
      }
      
     
      if(is.na(sum(B0[[1]]))!=1){
        B0_reduce=B0
        for(sizepointer in 1:d){
          B0_reduce=B0[[sizepointer]][1:p_reduce[sizepointer],1:r]
        }
      }

      if(is.na(sum(BB0[[1]]))!=1){
        BB0_reduce=resize_array(BB0,c(p_reduce))
      }
      if(is.na(sum(BB0_list[[1]][[1]]))!=1){
        Number_ini=length(BB0_list)
        BB0_list_reduce=BB0_list
        for(rep in 1:Number_ini){
          BB0_list_reduce[[rep]]=resize_array(BB0_list[[rep]],c(p_reduce))  
        }
        
      }    
      
    }else{
      #X_sample_reduce=X_sample
      r_reduce=r
    }  
  }


 test=list(X_sample_reduce=X_sample_reduce,r_reduce=r_reduce,B0_reduce=B0_reduce,B0_list_reduce=B0_list_reduce,BB0_list_reduce=BB0_list_reduce,BB0_reduce=BB0_reduce,shrink_factor=shrink_factor) 
return(test)
 }

warmstartupsize=function(warmstart=3,shrink_factor=1,p=1,X_sample_reduce_index=NA,BBreduce=NA){
  if(warmstart==1){
    if(shrink_factor<1){
      BBup=uparraysize(Xarray_new=BBreduce,pobjectup=p,index=X_sample_reduce_index[-length(X_sample_reduce_index)])
    }else{
      BBup=BBreduce  
    }
  }
  
  if(warmstart==3){
    if(shrink_factor<1){
      BBup=resize_array(BBreduce,p)  
    }else{
      BBup=BBreduce   
    }  
  }
  
  if(warmstart==4){
    BBup=upsizearray_blocker(BBreduce)  
  }
  
  if(warmstart==5){
    if(shrink_factor<1){
      BBup=resize_array(BBreduce,p)  
    }else{
      BBup=BBreduce 
    }
  }
  
  return(BBup)
}

warmstartprint=function(warmstart=3){
    if(warmstart==3){
      print("CP decomposition in the strategy of hua zhou's warmstart")  
    }
}

regionplot=function(bbeta){
  BB=full_R(bbeta)
  p=dim(BB)
  d=length(p)

  BB_plot=matrix(0,p[1],p[2])

  if(d==2){
    BB_plot=BB
  }
  if(d==3){
    for(i in 1:p[1]){
      for(j in 1:p[2]){
        BB_plot[i,j]=norm(as.matrix(BB[i,j,]),'f')
      }
    }
  }
  if(d==4){
    for(i in 1:p[1]){
      for(j in 1:p[2]){
        BB_plot[i,j]=norm(as.matrix(BB[i,j,,]),'f')
      }
    }
  }
  BB_plot=abs(BB_plot)
  BB_plot=BB_plot*1/max(BB_plot)
  return(BB_plot)
}



#return the plot by simple sum of square of the coefficients
regionplot_BB=function(BB,rough=0){
  p=dim(BB)
  d=length(p)
  if(d<=3){
    BB_plot=matrix(0,p[1],p[2])
  }
  if(d==4){
    BB_plot=array(0,c(p[1],p[2],p[3]))  
  }
  if(d==2){
    BB_plot=BB
  }
  if(d==3){
    for(i in 1:p[1]){
      for(j in 1:p[2]){
        BB_plot[i,j]=norm(as.matrix(BB[i,j,]),'f')
        if(rough==1){
          if(BB_plot[i,j]>0){
            BB_plot[i,j]=1
          }else{
            BB_plot[i,j]=0
          }
        }
      }
    }
  }
  if(d==4){
    for(i in 1:p[1]){
      for(j in 1:p[2]){
        for(k in 1:p[3]){
        BB_plot[i,j,k]=norm(as.matrix(BB[i,j,k,]),'f')
        if(rough==1){
        if(BB_plot[i,j,k]>0){
          BB_plot[i,j,k]=1
        }else{
          BB_plot[i,j,k]=0
        }
        }
        }
      }
    }
  }
  BB_plot=abs(BB_plot)
  BB_plot=BB_plot*1/max(BB_plot)
  return(BB_plot)
}

predictvalue_nonlinear=function(b0,BB,X,y,family='gaussian'){
  d=length(dim(X))-1
  if( family=='gaussian'){
    n=length(y)
    BB=as.tensor(BB)
    if(class(X)[1]=="Tensor"){
    }else{
      X=as.tensor(X)
    }
    yhat=c()
    for(i in 1:n){
      if(d==4){
        yhat[i]=innerProd(BB,X[,,,,i])+b0
      }
      if(d==3){
        yhat[i]=innerProd(BB,X[,,,i])+b0  
      }
    }
    rm(X)
    gc()
    
    
  }
  
  
  yy=c()
  yyhat=c()
  for(i in 1:n){
    if(y[i]>=0){
      yy[i]=1
    }else{
      yy[i]=-1
    }
    if(yhat[i]>=0){
      yyhat[i]=1
    }else{
      yyhat[i]=-1
    }
    
    
  }
  true=sum(yyhat==yy)
  return(true)
  
  }
#BB=BBbig_nonlinear[,,,,49]
#b0=b_validation_test_lambda_R_nonlinear[49,1]
#predictvalue_nonlinear(b0,BB,X,y,family='gaussian')
predictvalue_linear<-function(b0,BB,X,y,family='gaussian'){
  if( family=='gaussian'){
    n=length(y)
    if(class(BB)[1]=="Tensor"){
    }else{
      BB=as.tensor(BB)
    }
    
    
    if(class(X)[1]=="Tensor"){
    }else{
      X=as.tensor(X)
    }
    
    yhat=c()
    for(i in 1:n){
      yhat[i]=innerProd(BB,X[,,,i])+b0
    }
    yy=c()
    yyhat=c()
    for(i in 1:n){
      if(y[i]>=0){
        yy[i]=1
      }else{
        yy[i]=-1
      }
      if(yhat[i]>=0){
        yyhat[i]=1
      }else{
        yyhat[i]=-1
      }
      
      
    }
    true=sum(yyhat==yy)
    return(true)
  }
}
#predictvalue_linear(b0,BB,X_train,y_train,family='gaussian')





#to calculate  the l2 norm of the function
polyterm=function(alpha,beta,k,b){
if(beta==0){
  res=b^(alpha+1)/(alpha+1)
}else{
  res=b^(alpha+1)*(b+k)^beta/(alpha+1)-beta/(alpha+1)*polyterm(alpha+1,beta-1,k,b)
}
return(res)
}

polytermtrans=function(alpha,k1,beta,k2,b){
if(k2>=k1){
  alpha_middle=alpha
  k1_middle=k1
  alpha=beta
  k1=k2
  beta=alpha_middle
  k2=k1_middle
}
  k=k1-k2
  b=b-k1
  res=polyterm(alpha,beta,k,b)
  return(res)
}

#defaul is 6 basis
polyterm_one=function(knots,coef){
  K=length(coef)
  b=knots[K-1]
  a=knots[1]
  kksequence=c(0,0,0,knots[2:(K-2)])
  value_matrix=matrix(0,K,K)
  ordersequence=c(1,2,3,rep(3,(K-3)))
  for(i in 1:K){
    for(j in 1:K){
      if(i>=4||j>=4){
      value_matrix[i,j]=coef[i]*coef[j]*polytermtrans(ordersequence[i],kksequence[i],ordersequence[j],kksequence[j],b)
      }
    }
    
  }
  for(i in 1:3){
    for(j in 1:3){
      ordersmall=ordersequence[i]+ordersequence[j]+1
      value_matrix[i,j]=coef[i]*coef[j]*(b^(ordersmall)/ordersmall-a^(ordersmall)/ordersmall)
    }
  }
  
  
  return(sqrt(sum(value_matrix)))
}

#if using knots, it means  all the spline use the same knots. If using knotsBB, it means for 
#each  entry, the splien use different knots.
BB_l2norm_squre=function(BB,knots,knotsBB=NA,dimension=2){
  p=dim(BB)
  d=length(p)
  BB_plot=array(0,c(p[1:(d-1)]))
  if(d==4){
    for(i in 1:p[1]){
      for(j  in 1:p[2]){
        for(k in 1:p[3]){
          coef=BB[i,j,k,]
          if(is.na(knotsBB)==1){
          BB_plot[i,j,k]=polyterm_one(knots,coef)
          }else{
            BB_plot[i,j,k]=polyterm_one(knotsBB[i,j,k,],coef)  
          }
        }
      }
    }
  }
  
  if(d==3){
    for(i in 1:p[1]){
      for(j  in 1:p[2]){
        
          coef=BB[i,j,]
          if(is.na(knotsBB)==1){
          BB_plot[i,j]=polyterm_one(knots,coef)
        }else{
          BB_plot[i,j]=polyterm_one(knotsBB[i,j,],coef)
        }
        
      }
    }
  }
  
  if(dimension==2){
    if(d==4){
      #BB_plot_middle=BB_plot[,,1]/max(BB_plot[,,1])+BB_plot[,,2]/max(BB_plot[,,2])+BB_plot[,,3]/max(BB_plot[,,3])
      BB_plot_middle=BB_plot[,,1]^2+BB_plot[,,2]^2+BB_plot[,,3]^2
      
      BB_plot=sqrt(BB_plot_middle)
    }
  }
  
  rough=0
  if(rough==1){
    maxvalue=0.2*max(BB_plot)
    for(i in 1:p[1]){
      for(j in 1:p[2]){
        if(BB_plot[i,j]>maxvalue){
          BB_plot[i,j]=1
        }else{
          BB_plot[i,j]=0
        }
      }
    }
  }
  
  BB_plot=BB_plot/max(BB_plot)
  return(BB_plot)
  
}

BIC_computation<-function(){return(0)}
