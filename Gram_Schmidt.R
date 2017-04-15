## Function to perform Gram-schmidt factorization to obtain least square estimates
Gram=function(X,Y){
  ## X: Design matrix of training set
  ## Y: Responses vector of training set
  ## Output will be a list containing the least square estimates
  n=dim(X)[1]
  p=dim(X)[2]
  Q=matrix(rep(0,n*p),nrow=n,ncol=p)
  R=matrix(rep(0,p*p),nrow=p,ncol=p)
  v=X[,1]
  R[1,1]=sqrt(v%*%v)
  Q[,1]=v/R[1,1]
  for(j in 2:dim(X)[2]){
    v=X[,j]
    for(i in 1:(j-1)){
      R[i,j]=Q[,i]%*%X[,j]
      v=v-R[i,j]*Q[,i]
    }
    R[j,j]=sqrt(v%*%v)
    Q[,j]=v/R[j,j]
  }
  beta_hat=solve(R)%*%t(Q)%*%Y
  colnames(beta_hat)=c('Estimated coefficients')
  return(list(Estimates=beta_hat,UT=R,Orth=Q))
}
