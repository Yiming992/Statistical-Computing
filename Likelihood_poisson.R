NLPois<- function(Y,X,theta){
  X<-as.matrix(X)
  Y<-as.matrix(Y)
  Y<-as.vector(Y)
  
  theta2<- theta[2:length(theta)]
  alpha<-theta2[seq(1,length(theta2)-1,2)]
  beta<- theta2[seq(2,length(theta2),2)]
  j<- length(theta)
  J<- length(beta)
  likel<-NULL
  for (i in 1:(dim(X)[1])){
    if (X[i,1]==0){
      pi<-1}
    else{
      pi<-exp(alpha[X[i,1]]+X[i,2]*beta[X[i,1]])/(1+exp(alpha[X[i,1]]+X[i,2]*beta[X[i,1]]))
    }
    likel[i]<-pi
    likel<-as.vector(likel)
  }
  lnlklhd<- sum(-theta[1]*likel+Y*log(theta[1])+Y*log(likel)-log(factorial(Y)))
  library(phonTools)
  score<-zeros(j)
  score[1]<-sum(Y/theta[1]-likel)
  inform<- zeros(j,j)
  inform[1,1]<-sum(Y)/(theta[1]**2)
  for (m in 1:J){
    y<-Y[X[,1]==m]
    like<-likel[X[,1]==m]
    x<-X[X[,1]==m,2]
    
    score[2*m]<-sum((y-theta[1]*like)*(1-like))
    score[2*m+1]<- sum((y-theta[1]*like)*(1-like)*x)
    inform[1,2*m]<-sum(like*(1-like))
    inform[2*m,1]<-sum(like*(1-like))
    inform[1,2*m+1]<- sum(like*(1-like)*x)
    inform[2*m+1,1]<- sum(like*(1-like)*x)
    inform[2*m,2*m]<-sum((theta[1]-2*theta[1]*like+y)*like)
    inform[2*m,2*m+1]<-sum((theta[1]-2*theta[1]*like+y)*like*(1-like)*x)
    inform[2*m+1,2*m]<-sum((theta[1]-2*theta[1]*like+y)*like*(1-like)*x)
    inform[2*m+1,2*m+1]<- sum((theta[1]-2*theta[1]*like+y)*like*(1-like)*x^2)
  }
  return(list(lnlklhd=lnlklhd,score=score,inform=inform))
  
  
  
}

D<-read.table('http://www.ics.uci.edu/~dgillen/STAT211/Data/NLPoisData.txt')
NLPois(Y=D$Y,X=D[,-c(1,3)],theta=theta)

