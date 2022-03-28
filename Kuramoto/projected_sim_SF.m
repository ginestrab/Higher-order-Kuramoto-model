%close all;
clear all;
N = 200;
c=3;
Tmax = 10;
dt = 0.005; 
sigma_max=2;
omega = randn(N,1); %1.0*tan(pi*(rand(N,1)-0.5));% Lorenzian or Gaussian
omega = omega - mean(omega);%remove mean
Omega = 0.0; 



gamma = 2.5;
m=6;%power-law with min degree m, exp gamma
for i=1:N,
kas(i) = m*rand(1,1).^(-1/(gamma-1));
while kas(i)>sqrt(m*N),
    kas(i) = m*rand(1,1).^(-1/(gamma-1));
end
end


kave = mean(kas); %mean target degree
uno = ones(N,1); %vector of ones
clu = kas'*kas/(kas*uno); %matrix of probabilities
crand = rand(N,N);
a = crand<clu;
a=triu(a,1);
a=a+a';
k = sum(a)'; %degrees

omegaedge=Omega+rand(N,N);
omegaedge=triu(a.*omegaedge,1);
omegahat = sum(omegaedge,1)'-sum(omegaedge,2); %projected frequencies

thetas = 2*pi*rand(N,1); %initial conditions
psis = 2*pi*rand(N,1);

% increase sigma
j = 1; %counter
for sigma = 0:0.02:sigma_max;

    R0_ave = 0;
    R0hat_ave = 0;
    R1_ave = 0;
    %time evolution
    for t = 0:dt:Tmax
        r0 = sum(exp(1i*thetas))/N; %calculate R0,R_1^down, R_1
        r1 = sum(exp(1i*psis))/N;
        R0 = abs(r0);
        R1 = abs(r1);
        r0hat = sum(k.*exp(1i*thetas))/sum(k);
        R0hat = abs(r0hat);
   
        %evolve theta, psi using Eqs. (20), (24)
        thetap = omega +sigma*R1*imag(exp(-1i*thetas).*a*exp(1i*thetas));
        psisp = omegahat + sigma*R0*a*sin(psis) - sigma*R0*k.*sin(psis);
   
        thetas = thetas + dt*thetap; %forward Euler, works well in this case
        psis = psis + dt*psisp;
   
        if(t > (4*Tmax)/5) %time average in final fifth
            R0_ave = R0_ave + R0;
            R1_ave = R1_ave + R1;
            R0hat_ave = R0hat_ave + R0hat;
        end
   
    end
    sigma
    sigmasup(j) = sigma;
    ere1up(j) = R1_ave/(Tmax/5/dt);
    ere0up(j) = R0_ave/(Tmax/5/dt);
    ere0hatup(j) = R0hat_ave/(Tmax/5/dt);
    j = j+1;

    end

hold on
plot(sigmasup,ere1up,'o-') %plot Psi vs time
ylim([0,1])


j = 1;
for sigma = sigma_max:-0.02:0;

    R0_ave = 0;
    R0hat_ave = 0;
    R1_ave = 0;
    %time evolution
    for t = 0:dt:Tmax
        r0 = sum(exp(1i*thetas))/N; %calculate R0,R_1^down, R_1
        r1 = sum(exp(1i*psis))/N;
        R0 = abs(r0);
        R1 = abs(r1);
        r0hat = sum(k.*exp(1i*thetas))/sum(k);
        R0hat = abs(r0hat);
   
        %evolve theta, psi using Eqs. (20), (24)
        thetap = omega +sigma*R1*imag(exp(-1i*thetas).*a*exp(1i*thetas));
        psisp = omegahat + sigma*R0*a*sin(psis) - sigma*R0*k.*sin(psis);
   
        thetas = thetas + dt*thetap;
        psis = psis + dt*psisp;
   
        if(t > (4*Tmax)/5) %time average in final fifth
            R0_ave = R0_ave + R0;
            R1_ave = R1_ave + R1;
            R0hat_ave = R0hat_ave + R0hat;
        end
   
    end
    sigma
    sigmasdown(j) = sigma;
    ere1down(j) = R1_ave/(Tmax/5/dt);
    ere0down(j) = R0_ave/(Tmax/5/dt);
    ere0hatdown(j) = R0hat_ave/(Tmax/5/dt);
    j = j+1;


    end

%plot R_1^{down}
%figure
%plot(sigmasup,ere1up,'o-') %plot Psi vs time
%hold on
%plot(sigmasdown,ere1down,'o-') %plot Psi vs time
%ylim([0,1])

%plot R_0
%figure
%plot(sigmasup,ere0up,'o-') %plot Psi vs time
%hold on
%plot(sigmasdown,ere0down,'o-') %plot Psi vs time
%ylim([0,1])

%plot R_0hat
%figure
%plot(sigmasup,ere0hatup,'o-') %plot Psi vs time
%hold on
%plot(sigmasdown,ere0hatdown,'o-') %plot Psi vs time
%ylim([0,1])
