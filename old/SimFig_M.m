% this file sets the parameters and does the simulation for making the
% m-step fig.  stimulus is multidimensional gaussian noise. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data                  (smc_em_bern_real_exp)
% S:    simulation states for both samplers (smc_em_bern_main)
% M:    moments for both samplers           (smc_em_bern_main)
% E:    parameter estimates (not used)      (smc_em_bern_main)
% fig:  see fig file for details            (GetMFig1)

%% start function
clear; clc;

%% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.Nsec    = 20;               %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 1;                %number of spike history terms
Sim.freq    = 10;               %frequency of observations
Sim.StimDim = 5;                %stimulus dimensionality

%% make code prettier
Sim.K       = Sim.Nsec/Sim.dt;  %total # of steps
Sim.K_o     = Sim.K/Sim.freq;   %number of observations
Sim.tvec    = 0:Sim.dt:Sim.Nsec-Sim.dt; %time vector

%% set "real" parameters
P.k         = 2*1.6*ones(Sim.StimDim,1);%linear kernel
P.k(2:2:Sim.StimDim) = -1*P.k(2:2:Sim.StimDim);
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = 0.1;              %std of noise on calcium
P.sigma_o   = 20*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = 0;                %jump size for h after spike
P.tau_h     = 1;                %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%% do EM recursion
Sim.van     = false;            %not vanilla particle filtering
for runs=1:10
    % generate stimulus
    Sim.x       = rand(Sim.StimDim,Sim.K);
    Sim.x       = Sim.x.*repmat(sin(linspace(0,Sim.Nsec*pi,Sim.K)),Sim.StimDim,1);
    Sim.x(1,:)  = 1;

    % get "real" data
    R           = smc_em_bern_real_exp_v3(Sim,P);

    for iters=1:20
        [S M(runs)]     = smc_em_bern_FoBa_v1(Sim,R,P);
        E(runs,iters)   = smc_em_bern_mstep_v1(Sim,P,S,M);
        P               = E(runs,iters);
    end
    save('E','E')
end

%% make some figs
GetFig_M1(Sim,R,S,M,E)