# Function using EM algorithm to find the MLE for the mixture of point 
# mass at zero and a poisson(lambda)
EM=function(y,epsilon,iterations){
  # Input:
  # y: Observations vector
  # epsilon: convergence threshold
  # iterations: maximal number of iterations 
  
  # Output:
  # list containing estimated parameters of interests and number of iterations used
  
  old_theta=0.5
  old_lambda=0 # Set initial values
  n=length(y)
  j=0
  while(j<=iterations){
    p0=0
    p1=0
    p1_y=0
    # E-step
    for(i in 1:n){
      if(y[i]==0){
        p0=p0+old_theta/((old_theta)+(1-old_theta)
                         *dpois(y[i],old_lambda)) 
        p1=p1+(1-old_theta)*dpois(y[i],old_lambda)/
          ((old_theta)+(1-old_theta)*dpois(y[i],old_lambda)) 
        p1_y=p1_y+0
      }
      else{
        p0=p0+0
        p1=p1+1
        p1_y=p1_y+y[i]
      }
      # p0: Responsibilities of class zero
      # p1: Responsibilities of class one
      # p1_y: denominator used to update lambda in below M-step
    }
    # M-step update theta and lambda
    new_theta=p0/n
    new_lambda=p1_y/p1
    j=j+1
    # Test whether convergence is met
    if(sqrt((new_theta-old_theta)^2+(new_lambda-old_lambda)^2)<epsilon){break
    }
    old_theta=new_theta
    old_lambda=new_lambda
  }
  if( j> iterations ) cat( "Convergence not met before maximum iterations reached \n" )
  return(list(theta=new_theta,lambda=new_lambda,iter=j))
}