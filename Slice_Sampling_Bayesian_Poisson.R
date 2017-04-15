# Function to obtain Log-posterior in poisson regression, with normal(mu, tau) prior assumed for every regression parameter
getLogPost=function(x,y,currentBeta, mu, tau, allBeta){
  eta = x%*%allBeta;
  logLike = t(y)%*%eta - sum(exp(eta));
  logPrior= -((currentBeta-mu)^2) / ((2*tau^2));
  logPosterior = logLike + logPrior;
  return(logPosterior)
}

# assume normal prior with mean 0 and sd 10 for every regression parameter

w = 3; ##Step-size
m = 100;## maximum steps

# Function that performs one iteration of slice sampling
doSliceSampling <- function(x,y,currentBeta,j,mu,tau,allBeta){
  z = getLogPost(x,y,currentBeta, mu, tau, allBeta) - rexp(1, 1);
  # Stepping out to obtain a [L, R] range  
  u = runif(1);
  L = currentBeta - w*u;
  R = L + w;
  v = runif(1);
  J = floor(m*v);
  K = (m-1) - J;
  allBeta[j] = L;
  while(J>0 && z < (getLogPost(x,y,L,mu,tau,allBeta)) ){
    L = L - w;
    J = J - 1;
    allBeta[j] = L;
  }
  allBeta[j] = R;
  while(K>0 && z < (getLogPost(x,y,R,mu,tau, allBeta))){
    R = R+w;
    K = K-1;
    allBeta[j] = R;
  }
  # Shrinkage to obtain a new sample 
  u = runif(1);
  newBeta = L + u*(R-L);
  allBeta[j] = newBeta;
  while(z > (getLogPost(x,y,newBeta, mu, tau, allBeta))){
    if(newBeta < currentBeta){
      L = newBeta;
    }else{
      R = newBeta;
    }
    u = runif(1);
    newBeta = L + u*(R-L);
    allBeta[j] = newBeta;
  }
  return(newBeta) 
} 


#Function to implement above functions for slice sampling in poisson regression with given dataset.
Slice <- function(X,Y,nIterations,Burnin) {
  # X: design matrix
  # Y: response vector
  # nIterations: maximal number of slice sampling iterations
  # Burnin: burnin period of sampling
  # Output: a matrix containing all samples after the burnin period
  p=dim(X)[2]
  beta = matrix(rep(0, (nIterations+1)*(p)), ncol = (p), nrow=nIterations+1)	 
  for(i in 2:nIterations) {
    # allBeta is the current vector of beta's 
    allBeta = beta[i-1,]; 
    # cycle through beta's one at a time and update them keeping the other beta's fixed.
    for(j in 1:p) {
      # This is the current beta_j we want to update
      currentBeta = beta[i-1, j];
      # We pass the required information to the slice sampling function,
      # and obtain a new value beta_j at iteration "iter"
      beta[i, j] = doSliceSampling(X,Y,currentBeta, j, 0, 10, allBeta);
      # We also update allBeta since we have to use the most current
      # values when sampling the next beta.
      allBeta[j] = beta[i, j];
    }
  }
  return(beta[(Burnin+2):nIterations,])
}