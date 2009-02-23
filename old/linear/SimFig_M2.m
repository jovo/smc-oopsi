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
Sim.Nsec    = 60;              %# of sec
Sim.N       = 100;              %total number of particles
Sim.M       = 0;                %number of spike history terms
Sim.freq    = 10;               %frequency of observations
Sim.StimDim = 5;                %stimulus dimensionality
Sim.pf      = 3;


%% make code prettier
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem         = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                    %fix number of steps
end
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec-Sim.dt*rem;%time vector

%% set "real" parameters
P.k         = 8*1.6*ones(Sim.StimDim,1);%linear kernel
P.k(2:2:Sim.StimDim) = -1*P.k(2:2:Sim.StimDim);
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = 0.1;              %std of noise on calcium
P.sigma_o   = 20*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = -1;                %jump size for h after spike
P.tau_h     = 1;                %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%% do EM recursion
Sim.van     = false;            %not vanilla particle filtering

% generate stimulus
Sim.x       = 2*rand(Sim.StimDim,Sim.T);
% Sim.x(1,:)  = 0.5;

Nspikes     = Sim.T;%[100 500];% 1000 2000];
conv        = false;
Eold        = P;
Eold.k      = 0*Eold.k;
% Eold.omega  = 0*Eold.omega;
% Eold.tau_c  = 2*P.tau_c;
% Eold.beta   = 2*P.beta;
% Eold.sigma_c= 2*P.sigma_c;
% Eold.sigma_o= 2*P.sigma_o;

for secs=1:length(Nspikes)
    for runs=1:10

        % get "real" data
        R       = smc_em_bern_real_exp_v3(Sim,P);
        fprintf('# spikes= %d, rate=%g\n',sum(R.n), sum(R.n)/Sim.Nsec)
        spt     = find(R.n,Nspikes(secs));
        T(secs) = spt(end);
        Sim.Nsec= T(secs)*Sim.dt;
        
        % SMC-EM
%         while ~conv
            [S M(runs)]     = smc_em_bern_FoBaMoSuffStats_v2(Sim,R,Eold);
            E(secs,runs)   = smc_em_bern_mstep_v1(Sim,R,S,M,Eold);
            if norm(E(secs,runs).k./Eold.k)<1.01
                conv=true;
            end
            Eold = E(secs,runs);
            Eold.k
%         end
        save('Mstep','E','M')
    end
end

%% make some figs
GetFig_M3(Sim,R,S,M,E)