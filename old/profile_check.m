% this m-file is for the purposes of profiling
clear; clc;

[Sim P] = InitializeStuff;

%% set simulation parameters
Sim.Nsec    = 1;                       %# of sec
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;   %time vector
Sim.x       = ones(Sim.StimDim,Sim.T);  %stimulus
Sim.Mstep   = true;            %whether to estimate parameters
Sim.Rparams = true;            %whether to estimate rate governing parameters
Sim.Cparams = true;            %whether to estimate calcium transition parameters
Sim.SpikHis = true;            %whether to estimate the spike history parameters

R           = smc_em_bern_real_exp_v5(Sim,P);   %generate data
[S M]       = smc_em_bern_FoBaMo_v5(Sim,R,P);   %do forward-backward and get moments
Enew        = smc_em_bern_mstep_v2(Sim,S,M,P);  %maximize parameters