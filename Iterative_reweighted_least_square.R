## Funtion to perform iterative reweighted least square for logistic regression 

iterative=function(X,Y,epsilon,iter){
  ## X: Design matrix
  ## Y: Response matrix
  ## epsilon: stopping criterion
  ## iter: maximal number of iterations
  
  ## results of execution will be a list that contains the
  ## estimated betas and number of iterations used
  
  m=dim(X)[2]
  old=matrix(rep(0,m),nrow=m,ncol=1)
  i=0
  while(TRUE){
    w=diag(as.numeric(exp(X%*%old)/((1+exp(X%*%old))^2)))
    z=X%*%old+(Y*((1+exp(X%*%old))^2/exp(X%*%old))-(1+exp(X%*%old)))
    new=solve((t(X)%*%w%*%X))%*%t(X)%*%w%*%z
    i=i+1
    if (sqrt(t(new-old)%*%(new-old))<epsilon){
      break
    }
    if(i>iter){
      break
    }
    old=new
    
  }
  if( i > iter ) cat( "Convergence not met before maximum iterations reached \n" )
  return(list(beta_hats=new,iterations=i))
}