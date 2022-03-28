function [a,kappa,SC,count]=generate_SF_graph(N,gamma,m)
%%%%%Generate a SF graph%%%%%%%%
 % N number of nodes
%m minimum  degree
%gamma power-law exponent
%a adjacency matrix
for i=1:N,
    k(i)=m*rand(1)^(-1/(gamma-1));
    while(k(i)>sqrt(N))
        k(i)=m*rand(1)^(-1/(gamma-1));
    end
end
kappa=k;
x=rand(N,N);
x=x<((k'*k)/sum(k));
a=(tril(x,-1));
a=a+a';
      
%%%%%%%%Generates simplicial complex cell array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%SC{1} list of nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%SC{2} list of links %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%SC{3} list of triangles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
SC{1}=[1:N]';
[I,J,V]=find(tril(a));
count=0;
for in=1:numel(V),
    i=I(in);
    j=J(in);
SC{2}(in,1)=I(in);
SC{2}(in,2)=J(in);
for k=1:N,
    if((a(i,k)>0)&&(a(j,k)>0)&&(k>i))
        count=count+1;
        SC{3}(count,1)=j;
        SC{3}(count,2)=i;
        SC{3}(count,3)=k;
    end
end

end
