 
function [sigma,weight,mu] =  gaussian_mixture(data,k,epsilon, niterations)
% [gparams,memberships] =  gaussian_mixture(data,K,init_method,epsilon, 
%                                                niterations,plotflag, RSEED)
% INPUTS
%  data: N x d real-valued data matrix
%  K: number of clusters (mixture components)

%  epsilon: convergence threshold used to detect convergence
%  niterations (optional): maximum number of iterations to perform (default 500)

%  
%
% OUTPUTS
 
%           weight = N*k weight matrix
%           mu(k,:) = d-dimensional mean vector for kth component 
%           sigma{1,k} = d x d covariance vector for kth component
%  
%



% initialize with randomly selected parameters
[n,m]=size(data);
mu=data(randperm(n,k),:);
for j=1:k
    sigma{1,j}=cov(data);
end

iter=0;
alpha=ones(1,k)./k;
oldweight=zeros(n,k);
oldsigma=sigma;
oldmu=mu




for j=1:k
         density1(:,j)=alpha(1,j)*mvnpdf(data,mu(j,:),sigma{1,j});
         end
     
     oldloglike=sum(log(sum(density1,2)));




while iter<niterations
    
    
% perform E-step...
     for j=1:k
         pdf(:,j)=alpha(1,j)*mvnpdf(data,mu(j,:),sigma{1,j});
     end
     
     total= sum(pdf,2);
     for j=1:k
         weight(:,j)=pdf(:,j)./total; %memebership weight matrix
     end
     
         
    
% perform M-step...
     A=sum(weight,1);
     for j=1:k
         
         alpha(1,j)=A(1,j)/n;
         Ne(1,j)=A(1,j);%effective number of points to component k
         sigma{1,j}=zeros(m,m);
     end
     
     for j=1:k
         
         mu(j,:)=(weight(:,j).'*data)./Ne(1,j);
         for i=1:n
         dist=data(i,:)-mu(j,:);
         contri=weight(i,j)*(dist.')*dist;
         sigma{1,j}=sigma{1,j}+contri;
         end
         sigma{1,j}=sigma{1,j}./Ne(1,j);
     end
     
     

% compute log-likelihood and print to screen.....
         for j=1:k
         density(:,j)=alpha(1,j)*mvnpdf(data,mu(j,:),sigma{1,j});
         end
     iter=iter+1;
     loglike=sum(log(sum(density,2)));
     likel(1,iter)=loglike;
          

% check for convergence.....
 
 if loglike-oldloglike<epsilon
     iteration=[1:iter];
     plot(iteration,likel)
     xlabel('Iterations')
     ylabel('loglikelihood')
     title('change of loglikelihood with time')
         break
     else
         oldloglike=loglike;
         
 end
     
end


