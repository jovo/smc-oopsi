% this m-file generates the data and then plots the fig demonstrating how
% the two different sampling strategies differ.  It generates the following:
% 
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data                  (smc_em_bern_real_exp)
% S:    simulation states for both samplers (smc_em_bern_main)
% M:    moments for both samplers           (smc_em_bern_main)
% fig:  see fig file for details            (GetSamplFig1)

%% start function
clear; clc;

%% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.Nsec    = 2.01;%1.5;        %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 0;                %number of spike history terms
Sim.freq    = 5;                %relative frequency of observations

%% make code prettier
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem         = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                    %fix number of steps
end
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = 0:Sim.dt:Sim.Nsec-Sim.dt*(rem+1);%time vector

%% generate stimulus
Sim.x       = rand(1,Sim.T); 

%% set "real" parameters
P.k         = 1.8;%0.69;        %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = 0.1;              %std of noise on calcium
P.sigma_o   = 20*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = -3;               %jump size for h after spike
P.tau_h     = 0.8;              %decay rate for spike history terms
P.sigma_h   = 0.1;              %std of noise on spike history terms

%% get "real" data
R           = smc_em_bern_real_exp_v3(Sim,P);
sum(R.n), sum(R.n)/Sim.Nsec
% R.n         = zeros(1,Sim.T);   %spike times
% R.C         = zeros(1,Sim.T);   %initialize calcium
% epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
% spt         = [111 212 313 414 515];%[95 140];         %forced spike times
% R.n(spt)    = 1;                %force spikes
% 
% for t=2:Sim.T                   %update calcium
%     R.C(t)  = (1-Sim.dt/P.tau_c)*R.C(t-1) + P.beta*R.n(t) + epsilon_c(t);
% end
% R.O = R.C + P.sigma_o*randn(1,Sim.T);%add noise to observations

%% do EM recursion
for i=2:3
    i
    Sim.pf=i;
%     Sim.van=false;
%     if i==1                     %first do a backwards sampler
%         Sim.real=true;
%     else                        %then a vanilla sampler
%         Sim.real=false;
%     end
    [S{i} M(i)] = smc_em_bern_FoBaMo_v3(Sim,R,P);
    mse(i).n=sum(S{i}.w_b.*(S{i}.n-repmat(R.n,Sim.N,1)).^2);
    mse(i).C=sum(S{i}.w_b.*(S{i}.C-repmat(R.C,Sim.N,1)).^2);
%     mse(i).h=sum(S(i).w_b.*(S(i).h-repmat(R.C,Sim.N,1)).^2);
end
% [h pn]=ttest2(mse(1).n,mse(2).n);
% [h pC]=ttest2(mse(1).C,mse(2).C);
% % [h ph]=ttest2(mse(1).h,mse(2).h);
% [sum(mse(1).n) sum(mse(2).n) pn; ...
%     sum(mse(1).C) sum(mse(2).C) pC;...
%     ]%sum(mse(1).h) sum(mse(2).h) ph]
    


%% make some figs
GetFig_TestSamplers1(Sim,R,S,M)