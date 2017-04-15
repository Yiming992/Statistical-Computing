function[cluster,means,iter]=KMeans(data,k,r)
[n,m]=size(data);
cluster=cell(r,k);
means=zeros(1,r);
iter=zeros(1,r);
markers = {'+','o','*','.','x','s','d','^','v','>','<','p','h'};

for i=1:r
    % define initial centers
    ini=data(randperm(n,k),:);
    centers=zeros(k,m);
    iter(1,i)=0;
    Dist=cell(1,k);
    Means=[]
    while centers~=ini ;
        centers=ini;
        for j=1:k;
        Dist{1,j}=sqrt(sum((data-repmat(centers(j,:),n,1)).^2,2));
        
        end
        for l=1:k;                              
            cluster1=data(Dist{1,l}<=Dist{1,1},:);  
            for c=2:k;                                  
                                                        
                cluster2=data(Dist{1,l}<=Dist{1,c},:);
                cluster1=intersect(cluster1,cluster2,'rows');
            end
            cluster{i,l}=cluster1;
        end
        ini=zeros(k,m);
        Error=zeros(1,k);
        
        for h=1:k;
            [M,N]=size(cluster{i,h});
            ini(h,:)=mean(cluster{i,h},1);
            Error(1,h)=sum(sum((cluster{i,h}-repmat(ini(h,:),M,1)).^2,2))/M;
        end
        
        iter(1,i)=iter(1,i)+1;
        Means(1,iter(1,i))=sum(Error);
    end
    iteration=(1:iter(1,i));
    plot(iteration,Means,'Marker',markers{i})
    title('Mean-squared-error vs number of iterations') 
    hold on
    
  
    error=zeros(1,k);
    for e=1:k
        [o,p]=size(cluster{i,e});
        error(1,e)=sum(sum((cluster{i,e}-repmat(ini(e,:),o,1)).^2,2))/o;
    end
    means(1,i)=sum(error);
 end

        
        
       
        
        
                
                
            
        
            
        
        
        
        