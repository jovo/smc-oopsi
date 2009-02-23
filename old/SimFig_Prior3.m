% this m-file generates the data and then plots the fig demonstrating
% how stonger priors make the sampling more efficient
% observations. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data                  (smc_em_bern_real_exp)
% S:    simulation states for both samplers (smc_em_bern_FoBaMo)
% M:    moments for both samplers           (smc_em_bern_FoBaMo)
% fig:  see fig file for details

%% start function
clear; clc;

%% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.freq    = 10;               %relative frequency of observations
Sim.Nsec    = 0.75;              %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 1;                %number of spike history terms

%% make code prettier
Sim.T       = Sim.Nsec/Sim.dt;  %total # of steps
Sim.T_o     = Sim.T/Sim.freq;   %number of observations
Sim.tvec    = 0:Sim.dt:Sim.Nsec-Sim.dt;%time vector

%% set "real" parameters
P.k         = 1;                %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = .1;               %std of noise on calcium
P.sigma_o   = 20*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = -5;                %jump size for h after spike
P.tau_h     = 0.01;              %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%% get "real" data
% Sim.x       = ones(Sim.T,1);
% R           = smc_em_bern_real_exp_v3(Sim,P);
% sum(R.n)/Sim.Nsec/.35

R.n         = zeros(1,Sim.T);   %spike times
R.C         = zeros(1,Sim.T);   %initialize calcium
epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
spt         = [42 91 99];    %forced spike times
R.n(spt)    = 1;                %force spikes

for t=2:Sim.T                   %update calcium
    R.C(t)  = (1-Sim.dt/P.tau_c)*R.C(t-1) + P.beta*R.n(t) + epsilon_c(t);
end

Sim.van     = false;            %backwards sampler (ie, not vanilla)
R.O         = R.C + P.sigma_o*randn(1,Sim.T);%add noise to observations
priors      = [1 5 5];

%% do EM recursions
for h=1:2
    Sim.M=h-1;
    for n=1:3                       %for all 3 different priors
        Sim.x       = ones(Sim.T,1);%extize stimulus
        if n==2, Sim.x(spt(1):spt(end)) = 5*sin(pi:pi/(spt(end)-spt(1)):2*pi)+5;    %add appropriate pulese
        elseif n==3, Sim.x(spt)         = priors(n);
        end
        [S M(h,n)] = smc_em_bern_FoBaMo_v1(Sim,R,P);%do forward-backward and get moments
    end
end

%% make some figs
GetFig_Prior7(Sim,R,M,spt,priors)