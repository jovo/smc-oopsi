% this file sets the parameters and does the simulation for making the
% schematic fig.  stimulus is a sinusoid. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data      (smc_em_bern_real_exp)
% S:    simulation states       (smc_em_bern_FoBaMo)
% M:    moments                 (smc_em_bern_FoBaMo)
% fig:  see fig file for details(GetSchemFig1)

%% start function
clear; clc;

%% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.Nsec    = 0.65;             %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 0;                %number of spike history terms
Sim.freq    = 5;                %frequency of observations

%% make code prettier
Sim.T       = Sim.Nsec/Sim.dt;  %total # of steps
Sim.T_o     = Sim.T/Sim.freq;   %number of observations
Sim.tvec    = 0:Sim.dt:Sim.Nsec-Sim.dt; %time vector

%% generate stimulus
Sim.x       = 1+2*sin(linspace(-2*pi,2*pi,Sim.T));

%% set "real" parameters
P.k         = 1;                %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = .1;               %std of noise on calcium
P.sigma_o   = 1*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = 0;                %jump size for h after spike
P.tau_h     = 1;                %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%% get "real" data
R = smc_em_bern_real_exp_v3(Sim,P);

% R.n         = zeros(1,Sim.T);   %spike times
% R.C         = zeros(1,Sim.T);   %initialize calcium
% epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
% spt         = [17 72 89];       %forced spike times
% R.n(spt)    = 1;                %force spikes
% 
% for t=2:Sim.T                   %update calcium
%     R.C(t)  = (1-Sim.dt/P.tau_c)*R.C(t-1) + P.beta*R.n(t) + epsilon_c(t);
% end
% R.O = R.C + P.sigma_o*randn(1,Sim.T);%add noise to observations

%% do EM recursion
Sim.van     = false;            %not vanilla particle filtering
[S M]       = smc_em_bern_FoBaMo_v1(Sim,R,P);%do forward-backward and get moments

%% make some figs
GetFig_Schem3(Sim,R,S,M)