# This function uses the Gibbs sampler simulate from the posterior
# distribution of a normal model with unknown mean and variance. 
# The mean and variance are given normal(mu0,tau02) and Inv-chi2(nu0,sigma02) priors respectively. 

Gibbs = function(x,mu0,tau02,nu0,sigma02,nIterations) {
  n = length(x) # number of observations
  mu = 0 # starting state of mu
  sigma = 1 #starting state of sigma
  
  for(i in 2:nIterations) {
    # The conditional distribution of mu given sigma is normal(mu_n,
    # sigma_n) when we use a normal(mu0, tau02) prior for mu.
    mu_n = (mu0/tau02 + sum(x)/sigma[i-1]) / (1/tau02 + n/sigma[i-1]);
    sigma_n = sqrt(1/(1/tau02 + n/sigma[i-1]));
    mu[i] = rnorm(1, mean = mu_n, sd = sigma_n);
    # The conditional distribution of variance sigma2 is Inv-chi2(nu_n, sigma02_n)
    # We can first sample z form chi2(nu_n) distribution then use nu_n*sigma02_n/z
    # as a sample from the scaled Inv-Chi2 distribution.  
    nu_n = nu0+n;
    nu = (sum( (x - mu[i])^2 ))/n;
    sigma02_n = (nu0*sigma02+n*nu)/(nu0+n);
    z = rchisq(1, nu_n);
    sigma[i] = nu_n*sigma02_n/z;
  }
  return(list(mu = mu, sigma = sigma));
}