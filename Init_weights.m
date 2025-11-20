function [Popsize,weights,neighbour] = Init_weights(H,M,T)
[weights,Popsize] = Latin(H,M);
distance = zeros(Popsize,Popsize); 
neighbour=zeros(Popsize,T); 
for i=1:Popsize
        for j=i+1:Popsize
            A=weights(i,:);B=weights(j,:);
            distance(i,j)=(A-B)*(A-B)';  
            distance(j,i)=distance(i,j); 
        end
        [~,sindex]=sort(distance(i,:)); 
        neighbour(i,:)=sindex(1:T); 
end
end

function [W,N] = Latin(N,M)
    [~,W] = sort(rand(N,M),1);
    W = (rand(N,M)+W-1)/N;
end