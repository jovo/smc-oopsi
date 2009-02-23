% this file sets the parameters and does the simulation for making the
% schematic fig.  stimulus is a sinusoid. It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data     
% S:    simulation states      
% M:    moments                
% fig:  a schematic figure

%% start function
clear; clc;

[Sim P] = InitializeStuff;

% set figure specific parameters
Sim.Nsec    = 0.65;                     %# of sec
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem         = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                    %fix number of steps
end
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec-Sim.dt*rem;%time vector
Sim.N       = 200;

% generate stimulus
sin_init    = -2.2*pi;
Sim.x       = 1+2*sin(linspace(sin_init,sin_init+4*pi,Sim.T));

% get "real" data
R.n         = zeros(1,Sim.T);           %spike times
R.C         = P.C_init*ones(1,Sim.T);      %initialize calcium
epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
spt         = [24 89];                  %forced spike times
R.n(spt)    = 1;                        %force spikes

for t=2:Sim.T                           %update calcium
    R.C(t)  = (1-P.a)*R.C(t-1) + P.A*R.n(t) + P.a*P.C_0 + epsilon_c(t);
end
F_mu        = P.alpha*Hill_v1(P,R.C)+P.beta;        %compute E[F_t]
F_var       = P.gamma*Hill_v1(P,R.C)+P.zeta;    %compute V[F_t]
R.F         = F_mu+sqrt(F_var).*randn(1,Sim.T);%add noise to observations
R.F(R.F<0)  = eps;                      %observations must be non-negative

% do EM recursion
[S M]   = smc_em_bern_FoBaMo_v5(Sim,R,P);
fprintf('\n')

%% make a fig

%let O be only observations at sample times
O           = R.F.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
ONaNind     = find(~isfinite(O));
Oind        = find(isfinite(O));
O(ONaNind)  = [];

% do some x-axis stuff
xmin    = Oind(2);
xmax    = Oind(end-6);
xind    = xmin:xmax;
xs      = [Sim.tvec(xmin) Sim.tvec(xmax)]+Sim.dt;
spikeAX = [xs 0 1];
hidAX   = [xs 0 2];
xticks  = Sim.tvec(Oind(1))+Sim.dt:Sim.dt*Sim.freq*2:Sim.tvec(Oind(end));
xticks  = xticks-xticks(1);

% do some other stuff
cmin    = min(min(min(R.C(xind),min(M.Cbar(xind)-sqrt(M.Cvar(xind))))));
cmax    = max(max(max(R.C(xind),max(M.Cbar(xind)+sqrt(M.Cvar(xind))))));

gray    = [0.75 0.75 0.75];
col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8; ccol(ccol>1)=1;
ind     = Sim.T:-1:1;
sw      = 4.5;
fs      = 14;

figure(1), clf, Nsubs=4;
set(gcf, 'color', 'w');

% external stimulus
i=1; subplot(Nsubs,1,i), cla, hold on
plot(Sim.tvec,Sim.x,'k','LineWidth',2)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(Sim.x(xind)) max(Sim.x(xind))+1.])
ylab=ylabel({'External';'Stimulus'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% hidden states
i=i+1; subplot(Nsubs,1,i), cla, hold on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'Color',gray,'LineWidth',2)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis(hidAX)
ylab=ylabel({'True [Ca^2^+]';'and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',gray)

% observation state
i=i+1; subplot(Nsubs,1,i), cla, hold on
plot(Oind*Sim.dt,O,'.k','LineWidth',1,'markersize',7)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(O(2:end-6)) max(O(2:end-6))])
ylab=ylabel({'Observations'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% inference
i=i+1; subplot(Nsubs,1,i), cla, hold on,

% calcium
ptiles = GetPercentiles([.25 .75],S.w_b,S.C);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-cmin)/(cmax-cmin)+1,ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(M.Cbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')

% spikes
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
BarVar=M.nbar+M.nvar;
BarVar(BarVar>1)=1;
stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',sw)
stem(Sim.tvec,M.nbar,'Marker','none','Color',col(2,:),'LineWidth',sw)

% set stuff
set(gca,'YTick',[]), set(gca,'YTickLabel',[],'XTick',xticks),
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
axis(hidAX)
ylab=ylabel({'Inferred [Ca^2^+]';'and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(2,:))
xlab=xlabel('Time (sec)');

%% make eps
fig=figure(1);
wh=[7 3.5];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\schem