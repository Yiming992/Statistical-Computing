## Function to perform Householder factorization to obtain the least square estimates
House=function(X,Y){
  ## X: Design matrix of training set
  ## Y: Responses vector of training set
  ## Output will be a list containing the least square estimates
  n=dim(X)[1]
  p=dim(X)[2]
  U=matrix(rep(0,n*p),nrow=n,ncol=p)
  for(k in 1:p){
    w=X[k:n,k]
    w[1]=w[1]-sqrt(w%*%w)
    u=w/sqrt(w%*%w)
    U[k:n,k]=u
    X[k:n,k:p]=X[k:n,k:p]-2*u%*%(t(u)%*%X[k:n,k:p])
    Y[k:n]=Y[k:n]-2*u%*%(t(u)%*%Y[k:n])
  }
  R=X[1:p,1:p]
  beta_hat=solve(R)%*%Y[1:5]
  colnames(beta_hat)=c('Estimated coefficients')
  return(list(Estimates=beta_hat))
}