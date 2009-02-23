% this file sets the parameters and does the simulation for making the
% schematic fig.  stimulus is a sinusoid. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data      (smc_em_bern_real_exp)
% S:    simulation states       (smc_em_bern_FoBaMo)
% M:    moments                 (smc_em_bern_FoBaMo)
% fig:  see fig file for details(GetSchemFig1)

%% start function
clear; clc;

[Sim P] = InitializeStuff;

Sim.freq    = 5;
Sim.Nsec    = 10;               %# of sec
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem         = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                    %fix number of steps
end
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec-Sim.dt*rem;%time vector

Sim.M       = 1;
Sim.StimDim = 5;                %stimulus dimensionality
Sim.x       = rand(Sim.StimDim,Sim.T)*3-1;
% for t=1:Sim.freq:Sim.T
%     Sim.x(:,t:t+Sim.freq-1)=repmat(Sim.x(:,t),1,Sim.freq);
% end
Sim.x(1,:)  = ones(1,Sim.T);
P.k         = 1.3*ones(Sim.StimDim,1);%linear kernel
P.k(2:2:end)= P.k(2:2:end)-.5;
P.omega     = -6*P.k(1);
P.tau_h     = .025;                %decay rate for spike history terms
% P.gamma     = 1e-5;             %var gainh
% P.zeta      = 1e-5;            %var offset
P.C_0       = 1;
% P.A         = .1;
% P.tau_c     = .5;

R = smc_em_bern_real_exp_v5(Sim,P);

%%
%maximize tau_c, A, and C_0
A   = [R.C(1:end-1)*Sim.dt; -R.n(2:end); -repmat(Sim.dt,1,Sim.T-1)]';
b   = (R.C(2:end)-R.C(1:end-1))';
H = A'*A;
f = A'*b;
[minobsabc] = quadprog(H, f,[],[],[],[],[0 0 0],[inf inf inf]);
sigobs=sqrt(sum((R.C(2:end)-R.C(1:end-1)*(1-Sim.dt*minobsabc(1))-minobsabc(2)*R.n(2:end)-Sim.dt*minobsabc(3)).^2)/(Sim.T*Sim.dt));
[1/minobsabc(1), minobsabc(2), minobsabc(3)/minobsabc(1), sigobs; P.tau_c, P.A, P.C_0, P.sigma_c]
% Enew.tau_c  = 1/ve_x(1);
% Enew.A      = ve_x(2);
% Enew.C_0    = ve_x(3)/ve_x(1);

%% plot real n, C, F
figure(1), clf, hold on, title(num2str(sum(R.n)))
subplot(311), stem(Sim.tvec,R.n), ylabel('n')
subplot(312), plot(Sim.tvec,R.C), ylabel('[Ca++]'), axis([0 Sim.Nsec min(R.C) max(R.C)])
subplot(313), plot(Sim.tvec,R.F), ylabel('F_t'), axis([0 Sim.Nsec min(R.F) max(R.F)])
% subplot(312), plot(Sim.tvec,R.C-P.C_0), ylabel('[Ca++]'), axis([0 Sim.Nsec min(R.C-P.C_0) max(R.C-P.C_0)])
% subplot(313), plot(Sim.tvec,(R.F-P.beta)/P.alpha), ylabel('F_t'), axis([0 Sim.Nsec min((R.F-P.beta)/P.alpha) max((R.F-P.beta)/P.alpha)])

%%
figure(5),
subplot(511), plot(P.k'*Sim.x), ylabel('filt stim'), set(gca,'XTickLabel',[]), axis('tight'), title(['min=', num2str(min(P.k'*Sim.x)), ', max=',num2str(max(P.k'*Sim.x))])
subplot(512), plot(R.p), ylabel('p'), title(sum(R.n)/Sim.Nsec), set(gca,'XTickLabel',[])
subplot(513), stem(R.n), title(sum(R.n)), set(gca,'XTickLabel',[])
ISI=diff(find(R.n)); % x = 1:length(ISI); n=histc(ISI,x); n=n/sum(n); c=cumsum(n); bar(x,c),
subplot(514), plot(ISI), ylabel('ISI'), set(gca,'XTickLabel',[]), axis([0 length(ISI) 0 10]), title(length(find(diff(find(R.n))<Sim.freq))/sum(R.n))
subplot(515), plot(Sim.tvec,P.omega*R.h), ylabel('h'), axis('tight')

%% do EM recursion
Sim.pf      = 1;
Sim.Mstep   = true;%do M-step
E           = P;
Sim.Rparams = true;
E.k         = 0*E.k;
% Sim.Cparams = true;
E(1).tau_c  = 2*E(1).tau_c;
E(1).A      = 2*E(1).A;
E(1).C_0    = 2*E(1).C_0;
[S M]       = smc_em_bern_FoBaMo_v5(Sim,R,P);%do forward-backward and get moments
Enew        = smc_em_bern_mstep_v2(Sim,S,M,E);

% fprintf('\n')

%% % plot particles for n, C
% preset stuff for figs
gray    = [0.75 0.75 0.75];     %define gray
col     = [1 0 0; 0 .5 0];      %define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;%define colors for std
sw      = 2;                    %spike width
sh      = .2;                   %spike height
ind     = Sim.T:-1:1;           %inverse index for 'fill' plots

F           = R.F.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let F be only observations at sample times
FNaNind     = find(~isfinite(F));
Find        = find(isfinite(F));
F(FNaNind)  = [];
F(F<0)      = 0;
finv        = ((P.k_d.*(P.beta-F))./(F-P.beta-P.alpha)).^(1/P.n);

BarVar      = M.nbar+M.nvar; BarVar(BarVar>1)=1;

% plot particles
figure(2), clf,
subplot(121), hold on
stem(S.n'*sh,'Marker','none','color',ccol(2,:)),
stem(R.n*sh,'Marker','none','color',gray,'linewidth',2)
plot(sh*ones(size(Sim.tvec)),'k')
plot(repmat(1:Sim.T,Sim.N,1)',(S.C-P.C_0)'+2*sh,'color',ccol(2,:)),
plot(R.C-P.C_0+2*sh,'color',gray,'LineWidth',2)
plot(Sim.freq:Sim.freq:Sim.T,finv-P.C_0+2*sh,'.k','markersize',5)
plot(sum(S.w_f.*(S.C-P.C_0))+2*sh,'linewidth',2,'color',col(2,:))
axis([0 Sim.T 0 max(max(finv-P.C_0+2*sh),max(max(S.C-P.C_0))+2*sh)])
ylabel('[Ca++]')
title('particles')

% spike dist'n
subplot(122), hold on
stem(R.n*sh,'Marker','none','Color',gray,'LineWidth',sw)                  %plot true spikes
stem(BarVar*sh,'Marker','none','Color',ccol(2,:),'LineWidth',sw)          %plot forward spike var
stem(M.nbar*sh,'Marker','none','Color',col(2,:),'LineWidth',sw)        %plot forward spike mean
ptiles = GetPercentiles([.5 .95],S.w_b,S.C);
hfill=fill([1:Sim.T Sim.T:-1:1],[ptiles(1,:)-P.C_0 ptiles(2,ind)-P.C_0]+2*sh,ccol(2,:));
% hfill=fill([Sim.tvec Sim.tvec(ind)],[M.Cbar-sqrt(M.Cvar)-P.C_0 M.Cbar(ind)-P.C_0+sqrt(M.Cvar(ind))]+2*sh,ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(R.C-P.C_0+2*sh,'color',gray,'LineWidth',2)
plot(M.Cbar-P.C_0+2*sh,'linewidth',2,'color',col(2,:))
plot(Sim.freq:Sim.freq:Sim.T,finv-P.C_0+2*sh,'.k','markersize',5)
axis([0 Sim.T 0 max(max(finv-P.C_0+2*sh),max(max(S.C-P.C_0))+2*sh)])
title('mean +/- std')

%%
% figure(3), clf, imagesc(S.w_b), hold on, plot((1-max(S.w_b))*Sim.N,'k'), colorbar, title('forward-backward weights')
% figure(4), clf, imagesc(S.w_f), hold on, plot((1-max(S.w_f))*Sim.N,'k'), colorbar, title('forward-only weights')

