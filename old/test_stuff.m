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
Sim.Nsec    = 300;%1.5;              %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 1;                %number of spike history terms
Sim.freq    = 5;               %relative frequency of observations
Sim.StimDim = 1;                %dimensionality of stimulus

%% make code prettier
Sim.T       = Sim.Nsec/Sim.dt;  %total # of steps
Sim.T_o     = Sim.T/Sim.freq;   %number of observations
Sim.tvec    = 0:Sim.dt:Sim.Nsec-Sim.dt;%time vector

%% generate stimulus
Sim.x       = rand(Sim.T,1);     %make one input stationary

%% set "real" parameters
P.k         = 4.8;%0.69;             %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = 0.1;               %std of noise on calcium
P.sigma_o   = 20*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = [-0.1];% 0.05];                %jump size for h after spike
P.tau_h     = [0.8];%; 2];              %decay rate for spike history terms
P.sigma_h   = [0.1];% 0.1];             %std of noise on spike history terms

%% get "real" data
R           = smc_em_bern_real_exp_v3(Sim,P);
fprintf('T = %d, # of spikes = %d,  rate = %g\n', Sim.T, sum(R.n), sum(R.n)/Sim.Nsec)
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
Sim.van = true;
[S M]   = smc_em_bern_FoBaMo_v1(Sim,R,P);
% [S M]   = smc_em_bern_FoBaMoSuffStats_v2(Sim,R,P);
% Enew    = smc_em_bern_mstep_v1(Sim,R,S,M,P), P
fprintf('\n')
%% make some figs
% smc_em_bern_figs_v3(Sim,P,R,S,M)