function [Sim P] = InitializeStuff
% this function takes no input and outputs all the simulation parameters
% (embedded in the structure Sim.) and all the neuron parameters (embedded
% in the structure P.) necessary to simulate a neuron and do smc-em

% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.freq    = 5;                %observation frequnecy
Sim.N       = 100;              %total number of particles
Sim.M       = 0;                %number of spike history terms
Sim.pf      = 1;                %0: prior, 1: conditional
Sim.Nsec    = 10;               %# of sec
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   %time vector
Sim.StimDim = 1;                %# of stimulus dimensions
Sim.x       = ones(Sim.StimDim,Sim.T);  %stimulus

Sim.Mstep   = false;            %whether to estimate parameters
Sim.Rparams = false;            %whether to estimate rate governing parameters
Sim.Cparams = false;            %whether to estimate calcium transition parameters
Sim.SpikHis = false;            %whether to estimate the spike history parameters
Sim.conv    = false;            %whether the algorithm has converged yet
Sim.MaxIter = 25;               %max # of iters of EM

% set "real" parameters
P.k         = 1;                %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.A         = 0.1;              %jump size of calcium after spike
P.C_0       = 0;                %baseline calcium
P.C_init    = 0.05;             %initial calcium
P.sigma_c   = 0.01;             %std of noise on calcium

%spike history parameters
P.omega     = -1;               %jump size for h after spike
P.tau_h     = 1;                %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%fluorescence parameters
P.n         = 1.2;              %hill coefficient
P.k_d       = 1.3;              %dissociation constant
P.alpha     = 1;                %mean gain
P.beta      = 0;                %mean offset
P.gamma     = 1e-5;             %var gainh
P.zeta      = 1e-4;             %var offset

P.a         = Sim.dt/P.tau_c;   %for prettiness sake :-)