close all
%clear all
tic

%N = 100; %number of nodes
Tmax = 10; % max time to integrate
dt = 0.0003; %time step
sigma_max=5;
load('Results_Poisson_1000.mat')
incidence;
clear Rthetaup Rthetadown Xtheta Rpsiup Rpsidown
%load('data_NGFd3_N100_Gaussian.mat');
theta = 2*pi*rand(N,1); %initial conditions (random)

phi = 2*pi*rand(N2,1);


w = omega([1:N]); %Gaussian
w = w - mean(w);%enforce mean 0




what=omega([N+1:N+N2]);
what=what-sum(what);

%sweep sigma up
kcount = 1; %counter for k
B1t=transpose(B1);


for sigma = 0:0.03:sigma_max %sweep sigma up

    
    R1up_ave=0;
    Rthetaup_ave = 0;
    Rpsiup_ave = 0; %to store time averages
    Rpsi2up_ave = 0;

    te = 1; %time counter

    for t = 0:dt:Tmax %time loop


        psi=B1*phi;
        thetaB=B1t*theta;
        Xtheta = (exp(1i*thetaB));
        Xpsi = (exp(1i*psi));


        r0=abs(sum(exp(1i*theta)))/N;
        r1=abs(sum(exp(1i*phi)))/N2;
        r1down=abs(sum(Xpsi))/N;


        thetap = w - r1down*sigma*B1*imag(Xtheta);
        phip = what - sigma*r0*B1t*imag(Xpsi);
    
        theta = theta + dt*thetap; %forward Euler
        phi = phi + dt*phip;
    
    

    
        if(te > (4*Tmax/dt)/5) %time average in final fifth

            Rpsiup_ave = Rpsiup_ave + r1down;
            Rthetaup_ave = Rthetaup_ave + r0;
        end
        te = te +1;
    
    end

    R1up_ave = R1up_ave + r1;
    Rthetaup(kcount) = Rthetaup_ave/(te/5);
    Rpsiup(kcount) = Rpsiup_ave/(te/5);

    kup(kcount) = sigma; %store sigma
    kcount = kcount+1; %increase counter
    sigma
    end

    %plot

    plot(kup,Rthetaup,'+-');hold on; plot(kup,Rpsiup,'+-');%plot(kup,Rpsi2up,'+-');
    legend({'y = Rthetaup','y = Rpsiup'},'Location','northwest')
    ylim([0,1.1])
    hold on



    kcount = 1;
    for sigma = sigma_max:-0.03:0 %sweep sigma down, starting from previous final state
        R1down_ave=0;
        RpsiBdown_ave=0;
        Rthetadown_ave = 0;
        Rpsidown_ave = 0; %to store time averages

        te = 1;
        for t = 0:dt:Tmax
            psi=B1*phi;
            thetaB=B1t*theta;
            Xtheta = (exp(1i*thetaB));
            Xpsi = (exp(1i*psi));

            r0=abs(sum(exp(1i*theta)))/N;
            r1=abs(sum(exp(1i*phi)))/N2;
            r1down=abs(sum(Xpsi))/N;


            thetap = w - r1down*sigma*B1*imag(Xtheta); %Eqs (32) in notes
            phip = what - sigma*r0*B1t*imag(Xpsi); %- sigma*r1down*B2*imag(XpsiB);
    
            theta = theta + dt*thetap; %forward Euler
            phi = phi + dt*phip;
    

    
            if(te > 4*(Tmax/dt)/5) %time average in final fifth
                R1down_ave = R1down_ave + r1;
                Rpsidown_ave = Rpsidown_ave + r1down;
                Rthetadown_ave = Rthetadown_ave +r0;
            end
        te = te+1;
end

R1down(kcount)= R1down_ave/(te/5);
Rthetadown(kcount) = Rthetadown_ave/(te/5);
Rpsidown(kcount) = Rpsidown_ave/(te/5);
%RpsiBdown(kcount) = RpsiBdown_ave/(te/5);
kdown(kcount) = sigma; %store sigma
kcount = kcount+1; %increase counter
sigma
end

%plot
hold on
plot(kdown,Rthetadown,'+-');plot(kdown,Rpsidown,'+-');%plot(kdown,RpsiBdown,'+-')
legend({'y = Rthetadown','y = Rpsidown'},'Location','northwest')
ylim([0,1.1])

toc
