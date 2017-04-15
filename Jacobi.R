## Function to perform Jacobi algorithm to obtain the least square estimates
Jacobi=function(X,Y,k,epsilon){
  ## X: Design matrix of training set
  ## Y: Responses vector of training set
  ## k: Maximal number of iterations
  ## epsilon: Stopping criterion
  ## Output will be a list containing the least square estimates
  n=dim(X)[1]
  p=dim(X)[2]
  b=t(X)%*%Y
  A=t(X)%*%X
  P=diag(diag(A),nrow=p,ncol=p)
  I=diag(1,nrow=p,ncol=p)
  R_old=matrix(rep(0,p),nrow=p,ncol=1)
  i=1
  while(i<k){
    R_new=(I-solve(P)%*%A)%*%R_old+solve(P)%*%b
    error=sqrt(t(R_new-R_old)%*%(R_new-R_old))
    R_old=R_new
    i=i+1
    if(error<epsilon){
      break
    }
  }
  colnames(R_old)=c('Estimates')
  return(list(Estimates=R_old))
}