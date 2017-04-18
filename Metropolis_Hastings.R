# Assume normal prior with mean 0 and sd 100 for every regression parameters
# Function to obtain Log-posterior
getLogPost=function(x,y,currentBeta, mu, tau, allBeta){
eta = x%*%allBeta;
logLike = t(y)%*%eta - sum(exp(eta));
logPrior= -((currentBeta-mu)^2) / ((2*tau^2));
logPosterior = logLike + logPrior;
return(logPosterior)
}

# Function to implement metropolis-hastings in poisson regression
MH=function(X,Y,nIterations,Burnin){
  # X: design matrix
  # Y: response vector
  # nIterations: maximal number of MCMC iterations
  # Burnin: burnin period of MCMC sampling
  # Output: a matrix containing all samples after the burnin period

p=dim(X)[2]
beta = matrix(rep(0, (nIterations+1)*(p)), ncol = (p), nrow=nIterations+1);
for(i in 2:(nIterations+1)){
  beta_old=beta[i-1,]
  for(j in 1:p){#update regression parameter beta one at a time, while keeps the rest beta fixed
  beta_current=beta_old[j]# old beta
  old_post=getLogPost(X,Y,beta_current,0,100,beta_old)
  # Propose a new beta with normal(currentbeta,0.5)
  beta_new=rnorm(1,beta_current,0.5)
  beta_old[j]=beta_new
  new_post=getLogPost(X,Y,beta_new,0,100,beta_old)
  
  acceptance=min(1,exp((new_post+log(dnorm(beta_current,
                                           beta_new, 0.5)))-(old_post+log(dnorm(beta_new,beta_current,0.5)))))
  # Decide whether to accept the proposal
  u=runif(1)
  if(u>acceptance){
    beta_old[j]=beta_current
  }
  }
  beta[i,]=beta_old
}
  return(beta[(Burnin+2):nIterations,])
}
