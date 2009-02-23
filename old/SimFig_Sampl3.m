% this file sets the parameters and does the simulation for making the
% schematic fig.  stimulus is a sinusoid. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data      (smc_em_bern_real_exp)
% S:    simulation states       (smc_em_bern_main)
% M:    moments                 (smc_em_bern_main)
% fig:  see fig file for details(GetSchemFig1)

%% start function
clear; clc; [Sim P] = InitializeStuff;

%% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.Nsec    = 0.65;             %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 1;                %number of spike history terms
Sim.freq    = 10;               %frequency of observations

%% make code prettier
Sim.T       = Sim.Nsec/Sim.dt;  %total # of steps
Sim.T_o     = Sim.T/Sim.freq;   %number of observations
Sim.tvec    = 0:Sim.dt:Sim.Nsec-Sim.dt; %time vector

%% generate stimulus
Sim.x       = 1+sin(linspace(-2*pi,2*pi,Sim.T));

%% set "real" parameters
P.k         = 1;                %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = .1;               %std of noise on calcium
P.sigma_o   = 20*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = 0;                %jump size for h after spike
P.tau_h     = 1;                %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%% get "real" data
R           = smc_em_bern_real_exp_v3(Sim,P);

%% do EM recursion
for i=1:3
    Sim.i=i;
    [S{i} M(i)]       = smc_em_bern_FoBaMo_Test(Sim,R,P);%do forward-backward and get moments
    mse(i).n=sum(S{i}.w_b.*(S{i}.n-repmat(R.n,Sim.N,1)).^2);
    mse(i).C=sum(S{i}.w_b.*(S{i}.C-repmat(R.C,Sim.N,1)).^2);
end

%% make some figs
GetFig_Schem2(Sim,R,M)