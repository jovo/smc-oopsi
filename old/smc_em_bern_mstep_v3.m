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

[Sim P] = InitializeStuff;

% set simulation parameters
Sim.Nsec= 5;                       %# of sec
Sim.T   = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem     = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                %fix number of steps
end
Sim.T_o = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec= Sim.dt:Sim.dt:Sim.Nsec-Sim.dt*rem;%time vector

Sim.M       = 1;                    %# spike history terms
Sim.StimDim = 5;                    %stimulus dimensionality
Sim.x       = rand(Sim.StimDim,Sim.T);
% for t=1:Sim.freq:Sim.T              %make stimulus change at sampling rate
%     Sim.x(:,t:t+Sim.freq-1)=repmat(Sim.x(:,t),1,Sim.freq);
% end
Sim.x(1,:)  = ones(1,Sim.T);        %make a constant stim for baseline firing rate
P.k         = [1:Sim.StimDim]';     %linear kernel
P.k(2:2:end)= -1*P.k(2:2:end)+1;
P.omega     = -5*P.k(1);            %spike history weight
P.tau_h     = .025;                 %decay rate for spike history terms
P.gamma     = 1e-5;                 %var gainh
P.zeta      = 1e-4;                 %var offset
P.C_0       = 1;

R       = smc_em_bern_real_exp_v5(Sim,P);

%% plot real n, C, F
figure(1), clf, hold on, title(num2str(sum(R.n)))
subplot(311), stem(Sim.tvec,R.n), ylabel('n')
subplot(312), plot(Sim.tvec,R.C), ylabel('[Ca++]'), axis([0 Sim.Nsec min(R.C) max(R.C)])
subplot(313), plot(Sim.tvec,R.F), ylabel('F_t'), axis([0 Sim.Nsec min(R.F) max(R.F)])

%%
figure(5),
subplot(511), plot(P.k'*Sim.x), ylabel('filt stim'), set(gca,'XTickLabel',[]), axis('tight'), title(['min=', num2str(min(P.k'*Sim.x)), ', max=',num2str(max(P.k'*Sim.x))])
subplot(512), plot(R.p), ylabel('p'), title(sum(R.n)/Sim.Nsec), set(gca,'XTickLabel',[])
subplot(513), stem(R.n), title(sum(R.n)), set(gca,'XTickLabel',[])
ISI=diff(find(R.n)); % x = 1:length(ISI); n=histc(ISI,x); n=n/sum(n); c=cumsum(n); bar(x,c),
subplot(514), plot(ISI), ylabel('ISI'), set(gca,'XTickLabel',[]), axis([0 length(ISI) 0 10]), title(length(find(diff(find(R.n))<Sim.freq))/sum(R.n))
subplot(515), plot(Sim.tvec,P.omega*R.h), ylabel('h'), axis('tight')

Sim.Mstep   = true;                 %do M-step
Sim.conv    = false;                %not yet converged
Sim.i       = 1;                    %iteration #
Sim.MaxIter = 99;                   %max # of iters
E(1:Sim.MaxIter)= P;                %initialize parameter estimates

% Sim.Rparams     = true;
% E(1).k          = 0*E(1).k;
% E(1).omega      = 0;
% figure(8), clf, hold on, plot(P.k,'r')

Sim.Cparams     = true;
% E(1).tau_c      = 2*E(1).tau_c;
% E(1).A          = 2*E(1).A;
% E(1).C_0        = 2*E(1).C_0;

while Sim.conv==false;
    [S M]   = smc_em_bern_FoBaMo_v5(Sim,R,E(Sim.i));
    E(Sim.i+1) = smc_em_bern_mstep_v2(Sim,R,S,M,E(Sim.i));
    Sim.i   = Sim.i+1;
%     err     = norm(E(Sim.i).k-P.k);
%     fprintf('\nIteration #=%d, Norm of error=%g\n',Sim.i,err)
    dtheta  =  [(E(Sim.i).tau_c - E(Sim.i-1).tau_c)/E(Sim.i-1).tau_c;... 
                (E(Sim.i).A - E(Sim.i-1).A)/E(Sim.i-1).A;...
                (E(Sim.i).C_0 - E(Sim.i-1).C_0)/E(Sim.i-1).C_0];
    err  =  [(E(Sim.i).tau_c - P.tau_c)/P.tau_c;... 
                (E(Sim.i).A - P.A)/P.A;...
                (E(Sim.i).C_0 - P.C_0)/P.C_0];
    fprintf('\nIteration #=%d, percent estimate changes=%g, percent err=%g\n',Sim.i,dtheta,err)        
    save('M_stuff','Sim','R','P','S','M','E')
    plot(E(Sim.i).k), drawnow
    if Sim.i>=Sim.MaxIter           %check for convergence
        Sim.conv=true;
    end
end

% Nspikes     = Sim.T;%[100 500];% 1000 2000];
% conv        = false;
% Eold.omega  = 0*Eold.omega;
% Eold.tau_c  = 2*P.tau_c;
% Eold.beta   = 2*P.beta;
% Eold.sigma_c= 2*P.sigma_c;
% Eold.sigma_o= 2*P.sigma_o;

% for secs=1:length(Nspikes)
%     for runs=1:10
%
%         % get "real" data
%         R       = smc_em_bern_real_exp_v5(Sim,P);
%         fprintf('# spikes= %d, rate=%g\n',sum(R.n), sum(R.n)/Sim.Nsec)
%         spt     = find(R.n,Nspikes(secs));
%         T(secs) = spt(end);
%         Sim.Nsec= T(secs)*Sim.dt;
%
%         % SMC-EM
% %         while ~conv
%             [S M(runs)]     = smc_em_bern_FoBaMoSuffStats_v2(Sim,R,Eold);
%             E(secs,runs)   = smc_em_bern_mstep_v1(Sim,R,S,M,Eold);
%             if norm(E(secs,runs).k./Eold.k)<1.01
%                 conv=true;
%             end
%             Eold = E(secs,runs);
%             Eold.k
%         end
%         save('Mstep','E','M')
%     end
% end

%% make some figs
% O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
% Oind    = find(~isnan(O));
% xs      = [Sim.tvec(Oind(2)) Sim.tvec(Oind(end-1))];
% xmin    = find(Sim.tvec>xs(1),1);
% xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
% xind    = xmin:xmax;
%
% spikeAX = [xs 0 1];
% hidAX   = [xs 0 2];
%
% cmin=min(min(min(R.C(xind)),min(M.Cbar(xind)-sqrt(M.Cvar(xind)))),min(O(xind)));
% cmax=max(max(max(R.C(xind)),max(M.Cbar(xind)+sqrt(M.Cvar(xind)))),max(O(xind)));
%
% gray=[0.75 0.75 0.75];
%
% col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
% ccol=col+.8; ccol(ccol>1)=1;
% ind=Sim.T:-1:1;
%
% figure(5), clf, Nsubs=2;
% set(gcf, 'color', 'w');
%
% % linear kernel
% i=1; subplot(1,Nsubs,i), cla, hold on
% set(gca,'XTickLabel',[]),
% set(gca,'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(Oind)),
% set(gca,'YTick',[])
% axis([xs min(Sim.x(xind)) max(Sim.x(xind))+1.])
% ylab=ylabel({'Error'});
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
%
% % other parameters
% i=i+1; subplot(1,Nsubs,i), cla, hold on
% set(gca,'XTickLabel',[]),
% set(gca,'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(Oind)),
% set(gca,'YTick',[])
% axis(hidAX)
% ylab=ylabel({'Magnitude'});
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',gray)
%
% %%make eps
% fig=figure(6);
% wh=[6 3];
% set(fig,'PaperPosition',[0 11-wh(2) wh]);
% print -depsc C:\D\Research\liam\SMC_EM_GLM\Mstep