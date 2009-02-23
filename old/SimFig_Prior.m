% this m-file generates the data and then plots the fig demonstrating
% how stonger priors make the sampling more efficient 
% observations. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data                  (smc_em_bern_real_exp)
% S:    simulation states for both samplers (smc_em_bern_main)
% M:    moments for both samplers           (smc_em_bern_main)
% fig:  see fig file for details            (GetArrayFig1)

%% start function
clear; clc;

%% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.freq    = 10;               %relative frequency of observations
Sim.Nsec    = 0.8;              %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 1;                %number of spike history terms

%% make code prettier
Sim.K       = Sim.Nsec/Sim.dt;  %total # of steps
Sim.K_o     = Sim.K/Sim.freq;   %number of observations
Sim.tvec    = 0:Sim.dt:Sim.Nsec-Sim.dt;%time vector

%% set "real" parameters
P.k         = 1.0;              %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = .1;               %std of noise on calcium
P.sigma_o   = 10*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = 0;                %jump size for h after spike
P.tau_h     = 0.8;              %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%% get "real" data

R.n         = zeros(1,Sim.K);   %spike times
R.C         = zeros(1,Sim.K);   %initialize calcium
epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.K);%generate noise on calcium
spt         = [27 71 113];      %forced spike times
R.n(spt)    = 1;                %force spikes

for k=2:Sim.K                   %update calcium
    R.C(k)  = (1-Sim.dt/P.tau_c)*R.C(k-1) + P.beta*R.n(k) + epsilon_c(k);
end

Sim.van     = false;            %backwards sampler (ie, not vanilla)
R.O         = R.C + P.sigma_o*randn(1,Sim.K);%add noise to observations
priors      = [1.1062 2.3281 3.4814];%rate of [3 10 30] Hz

%% do EM recursions
for n=1:3
    Sim.x           = ones(Sim.K,1);
    Sim.x(spt)      = Sim.x(spt)*priors(n);
    [S(n) M(n)]     = smc_em_bern_FoBaMo_v4(Sim,R,P);
end

%% make some figs
GetFig_Prior2(Sim,R,S,M)