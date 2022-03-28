function [a,SC,count]=generate_random_graph(N,c)
%%%%%Generate a random graph%%%%%%%%
 % Number of nodes
%c average degree
%a adjacency matrix
x=rand(N,N);
x=x<(c/N);
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
