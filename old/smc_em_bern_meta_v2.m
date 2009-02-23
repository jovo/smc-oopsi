% this file is a meta-m-file, which calls various other m-files, depending
% on the desired use.
%
% inputs: this is not a proper function, so there are technically no inputs
% 
% outputs: again, this is not a proper function, noneless, the following structures are created:
% Sim:  simulation parameters
% R:    "real" neuron data
% P:    parameters of "real" neuron
% S:    simulation states
% E:    parameter estimates (from each EM iteration)
% M:    moments (from each EM iteration)

%% start function
clear; clc;
% dbstop if warning
warning off all

%% set simulation parameters
Sim.dt          = 0.005;                                %time step size (sec)
Sim.Nsec        = 0.8;                                  %# of sec
Sim.KernelSize  = 1;                                    %dimensionality of linear kernel for each neuron operating on stimulus
Sim.iters       = 1;                                    %number of EM iterations
Sim.N           = 100;                                  %total number of particles
Sim.M           = 1;                                    %number of spike history terms
Sim.dtSample    = 20*Sim.dt;                                %frequency of observations
Sim.pfs         = [1 2];                                %which particle filters to use: 1) prior, 2) one-step, 3) backwards

%% make code prettier
Sim.K           = Sim.Nsec/Sim.dt;                      %total # of steps
Sim.K_o         = Sim.Nsec/Sim.dtSample;                %number of observations
Sim.frac        = Sim.dtSample/Sim.dt;                  %# of samples per observation
Sim.tvec        = 0:Sim.dt:Sim.Nsec-Sim.dt;             %time vector

%% generate stimulus
Sim.x           = rand(Sim.KernelSize,Sim.K);           %the stimulus is white gaussian noise (we add the vector of ones to deal with the bias)
Sim.x(1,:)      = 1;                                    %make one input stationary
pulses          = linspace(10,12,4);
sp              = [55:56 75:76 95:96 115:116];
Sim.x(1,sp(1:2))= pulses(1);
Sim.x(1,sp(3:4))= pulses(2);
Sim.x(1,sp(5:6))= pulses(3);
Sim.x(1,sp(7:8))= pulses(4);

if Sim.KernelSize>1                                     %modulate time-varying stimuli with a sinusoid
        Sim.x(2:Sim.KernelSize,:) = Sim.x(2:Sim.KernelSize,:).*repmat(0.2*sin(linspace(-10*pi,10*pi,Sim.K)),Sim.KernelSize-1,1);%modulate input by a sinusoid
end

%% get "real" data
sumn=0;
while sumn~=2
    [R P] = smc_em_bern_real_exp_v1(Sim);
    sumn=sum(R.I(sp))
end

%% initialize EM parameters
E(1:max(Sim.pfs),1:Sim.iters+1) = P;                                       %initialize parameter estimates
% E.k = E.k+2*randn(Sim.KernelSize,1);        %add noise to initial conditions
% if Sim.M>0
%     E.omega = E.omega+E.omega.*rand(Sim.M,1);
% end

%% do EM recursion
clear S
for pf=Sim.pfs
    for iter=1:Sim.iters

        %% makes code prettier :-)
        B.sig2_c    = E(pf,iter).sigma_c^2*Sim.dt;
        B.sig2_o    = E(pf,iter).sigma_o^2;
        B.a         = 1-Sim.dt/E(pf,iter).tau_c;
        B.beta      = E(pf,iter).beta;
        B.kx         = E(pf,iter).k'*Sim.x;
        if Sim.M>0
            B.sig2_h    = E(pf,iter).sigma_h.^2*Sim.dt;
            B.g         = 1-Sim.dt/E(pf,iter).tau_h;
            B.omega     = E(pf,iter).omega;
        end

        %% forward step
        Sim.van = false;
        if pf==1
            S{pf} = smc_em_bern_FohBack_v4c(Sim,R,B);
%         elseif pf==2
%             S{pf} = smc_em_bern_FohBack_v4d(Sim,R,B);
        else
            Sim.van = true;
            S{pf} = smc_em_bern_FoVanil_v1(Sim,R,B);
        end

        %% backward step
        S{pf}.w_b = smc_em_bern_backwardPF_v4(Sim,S{pf},B);

        %% get moments and mse's
        M(pf,iter) = smc_em_bern_moments_v1(Sim,S{pf},R);

        %% maximization step
        E(pf,iter+1) = E(pf,iter);
        %     E(pf+1iter+1) = smc_em_bern_mstep_v1(Sim,S,E(,pf+1iter));
        %     norm(E(pf+1,iter+1).k-P.k)
    end
end

close all
GetSchemFig1G(Sim,R,S,B,M)
GetPriorFig2H(Sim,R,S)
