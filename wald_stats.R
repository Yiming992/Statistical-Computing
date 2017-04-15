mvWald<-function(model,cov.names){
  
  order<- match(cov.names,names(model$coefficients))
  
  estimates<-model$coefficients[order]
  estimates<-as.matrix(estimates)
  information<- solve(vcov(model)[order,order])
  wald_stats<- t(estimates)%*%information%*%estimates
  p_value<- pchisq(wald_stats,length(cov.names),lower.tail = FALSE)
  
  results<- cbind(wald_stats,p_value)
  rownames(results)<-c('results')
  colnames(results)<-c('wald_stats','p_value')
  return(results)
       
}


mvWald(glm1,c('wdata$k5One_under_six','wdata$k5two_or_more_under_six'))



linearContrast<- function(model,contr.coef,conf.level=0.95,exp.rslt=FALSE){
  est_contrast<- model$coefficients%*%t(t(contr.coef))
  sd<-sqrt(t(t(t(contr.coef)))%*%vcov(model)%*%t(t(contr.coef)))
  lower<- est_contrast-pnorm((1+conf.level)/2,0,1)*sd
  upper<- est_contrast+pnorm((1+conf.level)/2,0,1)*sd
  wald_stats<- (est_contrast/sd)^2
  p_value<- pchisq(wald_stats,1,lower.tail = FALSE)
  if(exp.rslt){
    lower<- exp(lower)
    upper<- exp(upper)
    est_contrast<-exp(est_contrast)
  }
  results<-cbind(est_contrast,sd,lower,upper,wald_stats,p_value)
  colnames(results)<-c('estimated_contrast','Standard_error','lower_bound','upper_bound','wald_statistics','p-value')
  row.names(results)<-c('results')
  
  return(results)
}