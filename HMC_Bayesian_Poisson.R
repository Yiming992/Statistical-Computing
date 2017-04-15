library(mvtnorm)

# Function that computes energy during hamiltonian iteration,with multivariate normal(mu, sigma) prior assumed for regression parameters
getU=function(x,y,mu,sigma,allBeta){
  allBeta=as.vector(allBeta)
  eta = x%*%allBeta;
  neg_logLike = -t(y)%*%eta +sum(exp(eta));
  neg_logPrior= (allBeta-mu)%*%solve(sigma)%*%(allBeta-mu)/2;
  U = neg_logLike + neg_logPrior;
  return(U)
}


# Function that implements Hamiltonian Monte Carlo in poisson regression with given dataset
Hamiltonian=function(X,Y,nIterations,Burnin,eps,L){
  # X: design matrix
  # Y: response vector
  # nIterations: maximal number of Hamiltonian iterations
  # Burnin: burnin period of sampling
  # eps:step size
  # L: max steps
  # Output: a matrix containing all samples after the burnin period
  m=dim(X)[2]
  beta = matrix(rep(0, (nIterations+1)*(m)), ncol = (m), nrow=nIterations+1)
  
  for(n in 2:(nIterations+1)){
    currentbeta=beta[n-1,]# current state
    current.p=as.vector(rmvnorm(1,sigma=diag(m)))#prpose a direction
    current.e=getU(X,Y,rep(0,m),diag(100,m),currentbeta)-dmvnorm(current.p, sigma=diag(m), log=TRUE)# current energy
    
    q=currentbeta
    p=current.p
    #update current state and direction
    p=p-eps/2*(-t(X)%*%Y+t(X)%*%exp(X%*%q)+1/2*(t(solve(diag(100,m)))+solve(diag(100,m)))%*%q)
    for(i in 1:L){
      q=q+eps*p
      if(i !=L){
        p=p-eps*(-t(X)%*%Y+t(X)%*%exp(X%*%q)+1/2*(t(solve(diag(100,m)))+solve(diag(100,m)))%*%q)
      }
    }
    p=p-eps/2*(-t(X)%*%Y+t(X)%*%exp(X%*%q)+1/2*(t(solve(diag(100,m)))+solve(diag(100,m)))%*%q)
    p=-p
    proposed.e=getU(X,Y,rep(0,m),diag(m),q)-dmvnorm(as.vector(p), sigma=diag(m), log=TRUE)# energy with updated state and direction
    acceptProb = min(1, exp(-proposed.e + current.e))
    # decide whether to accept new state
    if(is.finite(acceptProb) && runif(1)<acceptProb){
      beta[n,]=q 
    }else{
      beta[n,]=currentbeta
    }
  }
  return(beta[(Burnin+2):nIterations,])
}