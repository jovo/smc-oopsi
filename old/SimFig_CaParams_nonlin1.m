% this file sets the parameters and does the simulation for making the
% m-step fig.  stimulus is multidimensional gaussian noise. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data                  
% S:    simulation states for both samplers 
% M:    moments for both samplers           
% E:    parameter estimates (not used)      
% fig:  see fig file for details            

%% start function
clear; clc; close all

[Sim P] = InitializeStuff;

% set simulation parameters
Sim.Nsec= 10;                       %# of sec
Sim.T   = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem     = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                %fix number of steps
end
Sim.T_o = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec= Sim.dt:Sim.dt:Sim.Nsec-Sim.dt*rem;%time vector
Sim.M       = 1;                    %# spike history terms
Sim.StimDim = 5;                    %stimulus dimensionality
Sim.Mstep   = true;                 %do M-step
Sim.MaxIter = 99;                   %max # of iters
Sim.Ntrials = 10;
Sim.Cparams = true;
% Sim.Rparams = true;               %whether to estimate rate governing parameters

P.k         = .5*(1:Sim.StimDim)';     %linear kernel
P.k(2:2:end)= -1*P.k(2:2:end)+1;
P.omega     = -5*P.k(1);            %spike history weight
P.tau_h     = .025;                 %decay rate for spike history terms
P.tau_c     = 0.5;                  %decay rate of calcium
P.A         = 0.1;                  %jump size of calcium after spike
P.C_0       = 0.5;                  %baseline calcium
P.gamma     = 1e-5;                 %var gainh
P.zeta      = 1e-4;                 %var offset


% E(1:Sim.Ntrials)=EI;
% figure(2), clf, hold on, stem([P.tau_c; P.A; P.C_0; P.sigma_c],'r','linewidth',2), 

%%
for k=1:Sim.Ntrials
    fprintf('\n\nTrial #%d\n',k)
    Sim.x       = rand(Sim.StimDim,Sim.T);
    Sim.x(1,:)  = ones(1,Sim.T);        %make a constant stim for baseline firing rate
    R           = smc_em_bern_real_exp_v5(Sim,P);

    EI          = GetInitParams_v1(Sim,R.F,P);
    Enew        = EI;
    E(1)        = Enew;
    ErrNew      = norm([Enew.tau_c; Enew.A; Enew.C_0; Enew.sigma_c]-[P.tau_c; P.A; P.C_0; P.sigma_c]);
    ErrOld      = ErrNew;
    dtheta      = inf;
    i           = 0;                    %iteration number
    Sim.conv    = false;
    figure(1), clf, hold on, stem(([Enew.tau_c; Enew.A; Enew.C_0; Enew.sigma_c]-[P.tau_c; P.A; P.C_0; P.sigma_c])./([EI.tau_c; EI.A; EI.C_0; EI.sigma_c]-[P.tau_c; P.A; P.C_0; P.sigma_c]))
    while Sim.conv==false;
        [S M]   = smc_em_bern_FoBaMo_v5(Sim,R,Enew);
        Eold    = Enew;
        Enew    = smc_em_bern_mstep_v2(Sim,S,M,Enew);
        i       = i+1;
        ErrOld  = ErrNew;
        ErrNew  = norm([Enew.tau_c; Enew.A; Enew.C_0; Enew.sigma_c]-[P.tau_c; P.A; P.C_0; P.sigma_c]);
        dtheta  = norm([Enew.tau_c; Enew.A; Enew.C_0]-[Eold.tau_c; Eold.A; Eold.C_0])/norm([Eold.tau_c; Eold.A; Eold.C_0; Enew.sigma_c]);
        fprintf('\nIteration #%d, Norm of error=%g, derr=%g, dtheta=%g\n',i,ErrNew,ErrOld-ErrNew,dtheta)
        save('Ca_params','Sim','R','P','M','E')
        figure(1), stem(([Enew.tau_c; Enew.A; Enew.C_0; Enew.sigma_c]-[P.tau_c; P.A; P.C_0; P.sigma_c])./([EI.tau_c; EI.A; EI.C_0; EI.sigma_c]-[P.tau_c; P.A; P.C_0; P.sigma_c])), drawnow
        if i>Sim.MaxIter || dtheta<.006
            Sim.conv=true;
            if ErrNew>ErrOld
                E(k)=Eold;
            else
                E(k)=Enew;
            end
%             figure(2), plot([E(k).tau_c; E(k).A; E(k).C_0; E(k).sigma_c],'b','linewidth',2), drawnow
        end
%         figure(5), clf, hold on, plot(S.C'), plot(R.C,'r','linewidth',2), drawnow
    end
end
%%
% plot real n, C, F
figure(3), clf, hold on, title(num2str(sum(R.n)))
subplot(311), plot(Sim.tvec,R.p/Sim.dt), ylabel('n')
subplot(312), plot(Sim.tvec,R.C), ylabel('[Ca++]'), axis([0 Sim.Nsec min(R.C) max(R.C)])
subplot(313), plot(Sim.tvec,R.F), ylabel('F_t'), axis([0 Sim.Nsec min(R.F) max(R.F)])
% subplot(312), plot(Sim.tvec,R.C-P.C_0), ylabel('[Ca++]'), axis([0 Sim.Nsec min(R.C-P.C_0) max(R.C-P.C_0)])
% subplot(313), plot(Sim.tvec,(R.F-P.beta)/P.alpha), ylabel('F_t'), axis([0 Sim.Nsec min((R.F-P.beta)/P.alpha) max((R.F-P.beta)/P.alpha)])

%
figure(4),
subplot(511), plot(P.k'*Sim.x), ylabel('filt stim'), set(gca,'XTickLabel',[]), axis('tight'), title(['min=', num2str(min(P.k'*Sim.x)), ', max=',num2str(max(P.k'*Sim.x))])
subplot(512), plot(R.p), ylabel('p'), title(sum(R.n)/Sim.Nsec), set(gca,'XTickLabel',[])
subplot(513), stem(R.n), title(sum(R.n)), set(gca,'XTickLabel',[])
ISI=diff(find(R.n)); % x = 1:length(ISI); n=histc(ISI,x); n=n/sum(n); c=cumsum(n); bar(x,c),
subplot(514), plot(ISI), ylabel('ISI'), set(gca,'XTickLabel',[]), axis([0 length(ISI) 0 10]), title(length(find(diff(find(R.n))<Sim.freq))/sum(R.n))
subplot(515), plot(Sim.tvec,P.omega*R.h), ylabel('h'), axis('tight')

figure(5), clf, hold on,
plot(S.C')
plot(R.C,'r','linewidth',2)