
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%/*This program sets the parameters of used in the simulations*/ %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
s=-1; %flavor s=-1,0,1
d=1; % do not change
c=2.5; %Average degree for Poisson network
m=1; %Minimum degree of the power-law
gamma=2.5; % power-law exponent for power-law network%
N=100; % Total number of nodes%
Teq=10*N; %Equilibration time%
Ttot=100;
T=Ttot*N; %Total time%
omega0=2;% non zero frequency for the links%

t = linspace(0,Ttot,T);

%%%%%%Set model choice between random graph and scale-free network%%%%%%%%%%%
%count=0;
%while (count==0)
%[a,SC,count]=generate_random_graph(N,c) %Poisson network
%returns the adjacency matrix a and the simplicial complex cell array SC to be read by Boundary%
%%and the simplicial complex cell array SC to be read by Boundary
%[a,kappa,SC,count]=generate_SF_graph(N,gamma,m) %Uncorrelated SF network
%returns the adjacency matrix a and the simplicial complex cell array SC to be read by Boundary%
%%and the simplicial complex cell array SC to be read by Boundary

%end
    
%%%%%%%%%%%%% /*Set models for SC */%%%%%%%%%%%%%%%
%%%%%NGF%%%%%%%%%%%%%%
[a,r,kn,kt,SC] = NGF_d3(N,s,0,0);
%%%%%%%%%Configuration model of simplicial complex
%! cc ensembles_simpliciald3.c -lm
%! ./a.out
%clear SC
%SC{1}=load('SC0.txt');
%SC{2}=load('SC1.txt');
%SC{3}=load('SC2.txt');
%SC{4}=load('SC3.txt');
%b=load('SC1.txt');
%b(:,3)=ones(1,numel(b(:,1)));
%a=spconvert(b);
%a(N,N)=0;

%%%%%%%%%%%%%Prepare boundary matrixes and parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
a=(a+a')>0;
%check a is simple and symmetric%
k=sum(a,2); %degree sequence%
kav=sum(k)/N; %average degree%

B2=Boundary(SC{d+2},SC{d+1});
N2=numel(sum(B2,2));
N3=numel(sum(B2,1));
if(d>0)
B1=Boundary(SC{d+1},SC{d});
end


%%%%%%%%% Set the proper frequencies of node omega([1:N]) and the proper frequencies of links omega([N+1:N+N2]) Gaussian distributed the nodes have average frequency zero the links omega0
%%%%%%%%%%%%%%%%
   
l=1:(N+N2);
omega1=omega0-omega0*(l<(N+1))';
%omega1=omega0;
omega=omega1+randn(N+N2,1);
   
%%%%%%%%%%Save data in a mat file%%%%%%%%%%%%%%%%%
save('data3_NGF_N100_conf_Gaussian.mat');

  
