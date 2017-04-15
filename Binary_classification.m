function [weights] = logistic_train(data,labels,epsilon,maxiterations,SGflag,M)
%% Function to implement numerical optimization in logistic regression setting
%% data: training data set
%% labels: class labels in training set
%% epsilon: Convergence condition
%% maxiterations: maximum number of iterations
%% SGFflag: Indicator, if =1, performs gradient descent, otherwise performs Newton-raphson
%% M: bstch size used during the optimization
%% 
 [n,m] = size(data);
 data=[data ones(n,1)];
 [n,m]=size(data);
 weights= zeros(m,1);
 i=0;
 dist=1;
 accuracy=[];
 log_loss=[];
 lambda=0.001;
 if SGflag==1
     tic;
     while i<maxiterations && dist> epsilon;
     
     order= randperm(n,M);
     Data= data(order,:);
     Labels= labels(order,:);
     p= exp(data*weights)./(1+exp(data*weights));
     q= exp(Data*weights)./(1+exp(Data*weights));
     gradient= transpose(Data)*(Labels-q);
     weights= weights+lambda*gradient;
     prob= exp(data*weights)./(1+exp(data*weights));
     dist= sum(abs(prob-p))/n;
     i=i+1;
     chat=prob>0.5;
     accuracy(1,i)=100*sum(labels==chat)/length(labels);
     
     end
     i

     
 else
     
     while i<maxiterations && dist>epsilon;
 p= exp(data*weights)./(1+exp(data*weights));
 
 v= 1-p;
 di= p.*v;
 w= diag(di);
 hessian= -(transpose(data)*w*data)+eye(m)*10^(-6);
 gradient= transpose(data)*(labels-p);
 weights= weights-inv(hessian)*gradient;
 prob= exp(data*weights)./(1+exp(data*weights));
 dist= sum(abs(prob-p))/n;
 i=i+1;
 chat=prob>0.5;
 accuracy(1,i)=100*sum(labels==chat)/length(labels);
 log_loss(1,i)= -(transpose(labels)*log(prob)+transpose(1-labels)*log(1-prob));
     end
 i
 end

 
end
