function B = Boundary( SimplicesN, SimplicesN_1 )
%BOUNDARY Simplicial chain boundary matrix
%   Each row of SimplicesN represents an N-simplex
%   Each row of SimplicesN_1 represents an (N-1)-simplex
%   Each simplex is represented by a list of its vertices
%   Vertices are integers 1 to n
%
%   B(i,j)=1 (resp. -1) if the ith (N-1)-simplex is incident with the 
%   jth N-simplex with positive (resp. negative) orientation
%   
%   If N=1 or N=0 then B=[] (zero boundary)

%sort rows (i.e. assume standard orientation of simplices)
SimplicesN = sort(SimplicesN,2);
SimplicesN_1 = sort(SimplicesN_1,2);

%dimension, number of simplices, initialize B
DIM=size(SimplicesN,2); 
if DIM<2
    B=[];
    return;
elseif DIM-1~=size(SimplicesN_1,2) 
    error('ERROR: Simplex dimensions don''t match');
end
N1=size(SimplicesN,1);   %number of N-simplices
N2=size(SimplicesN_1,1); %number of (N-1)-simplices
B=zeros(N2,N1); 

%remove one column at a time, find N-1 simplex, assign sign
for i=1:DIM
    Baux = SimplicesN(:,setdiff(1:DIM,i));
    [~,ind] = ismember(Baux, SimplicesN_1, 'rows');
    if ind==0
        error('Boundary simplex not found.');
    end
    w=(-1)^(i-1);   %weight
    for j=1:N1
        B(ind(j),j)=w;
    end
end
    
end

