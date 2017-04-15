## Function to perform Newton method for logistic regression 

Newton=function(X,Y,lambda,beta,alpha,epsilon,iter){
  ## X: Design matrix
  ## Y: Response matrix
  ## lambda: initial step size
  ## beta: Parameter for backtracking, between 0 and 1
  ## alpha: Parameter for backtracking, between 0 and 0.5
  ## epsilon: stopping criterion
  ## iter: maximal number of iterations
  ## objective of this function is to minimize the negative log-likelihood
  
  ## results of execution will be a list that contains the
  ## estimated betas and number of iterations used
  
  m=dim(X)[2]
  old=matrix(rep(0,m),ncol=1,nrow=m)
  dist=1
  i=0
  while(i<=iter){
    old_likelihood=sum(log(exp(X%*%old)+1))-t(Y)%*%X%*%old
    p=exp(X%*%old)/(1+exp(X%*%old))
    w=diag(as.numeric(p*(1-p)))
    hessian=(t(X)%*%w%*%X)
    gradient=t(X)%*%(p-Y)
    t=lambda
    while(TRUE){## backtracking to select step-size
      update=old-t*solve(hessian)%*%gradient
      new_likelihood=sum(log(exp(X%*%update)+1))-t(Y)%*%X%*%update
      dist=old_likelihood-new_likelihood
      if(dist>=alpha*t*t(gradient)%*%solve(hessian)%*%gradient){
        break
      }
      t=beta*t
    }
    i=i+1
    old=update
    if(dist<epsilon){
      break
    }
  }
  if( i > iter ) cat( "Convergence not met before maximum iterations reached \n" )
  return(list(Beta_hats=update,iterations=i))
}
