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
Sim.M       = 1;

% generate stimulus
sin_init    = -2.2*pi;
Sim.StimDim = 5;
Sim.x       = .02*rand(Sim.StimDim,Sim.T);
Sim.x       = Sim.x.*(25*repmat(sin(linspace(sin_init,sin_init+4*pi,Sim.T)),Sim.StimDim,1));

% linear filter
P.k         = ones(Sim.StimDim,1);%linear kernel
P.omega     = -3;               %jump size for h after spike
P.tau_h     = .01;                %decay rate for spike history terms

% get "real" data
kx        = P.k'*Sim.x;                                 %external input to neuron
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
U_sampl   = rand(1,Sim.T);                              %generate random number to use for sampling
R.C       = P.C_init*ones(1,Sim.T);                              %initialize calcium
R.p         = zeros(1,Sim.T);                       %extize p_t because it must be updated iteratively
R.n         = zeros(1,Sim.T);                       %extize n_t because it must be updated iteratively
spt         = [24 89];                  %forced spike times
R.n(spt)    = 1;                        %force spikes
R.h         = zeros(Sim.M,Sim.T);                   %extize spike history because it must be updated iteratively
R.y         = zeros(1,Sim.T);
epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(Sim.M,Sim.T); %generate noise on spike history
for t=2:Sim.T                                       %update states
    R.h(:,t)= (1-Sim.dt./P.tau_h).*R.h(:,t-1)+R.n(t-1) + epsilon_h(:,t);%update h terms
    R.y(t)=kx(t)+P.omega'*R.h(:,t);                    %generate operand for rate function
    R.p(t)=1-exp(-exp(R.y(t))*Sim.dt);                 %generate rate
end %time loop

epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
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
M.pbar = sum(S.w_b.*S.p,1);
M.pvar = sum((repmat(M.pbar,Sim.N,1)-S.p).^2)/Sim.N;
pmin    = min(min(min(R.p(xind),min(M.pbar(xind)-sqrt(M.pvar(xind))))));
pmax    = max(max(max(R.p(xind),max(M.pbar(xind)+sqrt(M.pvar(xind))))));

M.hbar = sum(S.w_b.*S.h,1);
M.hvar = sum((repmat(M.hbar,Sim.N,1)-S.h).^2)/Sim.N;
hmin    = min(min(min(P.omega*R.h(xind),P.omega*min(M.hbar(xind)-sqrt(M.hvar(xind))))));
hmax    = max(max(max(P.omega*R.h(xind),P.omega*max(M.hbar(xind)+sqrt(M.hvar(xind))))));

cmin    = min(min(min(R.C(xind),min(M.Cbar(xind)-sqrt(M.Cvar(xind))))));
cmax    = max(max(max(R.C(xind),max(M.Cbar(xind)+sqrt(M.Cvar(xind))))));

gray    = [0.75 0.75 0.75];
col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8; ccol(ccol>1)=1;
ind     = Sim.T:-1:1;
sw      = 4;
lw      = 2;
fs      = 14;
tfs     = 18;

figure(1), clf, 
Nrows=7;
Ncols=2;
set(gcf, 'color', 'w');

% external stimulus
subplot('Position',[0.35 0.88 0.336 0.1]), cla, hold on
imagesc(Sim.x)
colormap('gray')
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xmin xmax 0.5 5.5])
ylab=ylabel({'External'; 'Stimulus'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% filtered stimulus
subplot('Position',[0.35 0.75 0.336 0.1]), cla, hold on, box on
plot(Sim.tvec,P.k'*Sim.x,'color','k','linewidth',lw)
set(gca,'XTick',xticks,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
set(gca,'YTick',[])
% set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(P.k'*Sim.x) max(P.k'*Sim.x)])
ylab=ylabel({'Linearly Filtered'; 'Stimulus'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% ylab=ylabel({'External';'Stimulus'});
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% title('True','fontsize',tfs);

% y_t
% ii=3; i=ii; subplot(Nrows,Ncols,i), cla, hold on, box on
% plot(Sim.tvec,R.y,'Color',gray,'LineWidth',lw)
% axis([xs min(R.y) max(R.y)])
% set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])

% P(spike) 
ii=5; i=ii; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,(R.p-pmin)/(pmax-pmin),'Color',gray,'LineWidth',lw)
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
% title('$p_t=e^{f(b+k^T x_t+\omega h_{t})dt}$','Interpreter','latex');

% spike train
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
% title('$n_t \sim p_t$','Interpreter','latex');

% spike history
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,(P.omega*R.h-hmin)/(hmax-hmin),'Color',gray,'LineWidth',lw)%(R.h-hmin)/(hmax-hmin)
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
% title('$h_{t} - h_{t-1} = - \frac{dt}{\tau_{h}} h_{t-1} + n_{t-1} + \sigma_{h} \sqrt{dt} \varepsilon_{t}$','Interpreter','latex');

% calcium
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray,'LineWidth',lw)%(R.C-cmin)/(cmax-cmin)
% set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs 0 1])
% title('$[\textrm{Ca}]_t - [\textrm{Ca}^{2+}]_{t-1} = - \frac{dt}{\tau_c} ([\textrm{Ca}^{2+}]_{t-1}-[\textrm{Ca}^{2+}]_0) + A n_t + \sigma_c \sqrt{dt} \varepsilon_{t}$','Interpreter','latex');
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})

% observation state
i=i+2.5; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Oind*Sim.dt,O,'.k','LineWidth',1,'markersize',7)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(O(2:end-6)) max(O(2:end-6))])
xlabel('Fluorescence Observations');
% ylab=ylabel({'Fluorescence'; 'Observations'});
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% title('$F_t = \alpha  S([\textrm{Ca}^{2+}]_t) + \beta + \sigma_{Ft} \epsilon_t$','Interpreter','latex');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y_t
% i=ii+1; subplot(Nrows,Ncols,i), cla, hold on, box on
% plot(Sim.tvec,R.y,'Color',col(2,:),'LineWidth',lw)
% set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
% ylab=ylabel({'Input to      ';'Neuron      '});
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% axis([xs min(R.y) max(R.y)])

% P(spike)
i=ii+1; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,S.p);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-pmin)/(pmax-pmin),ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(R.p-pmin)/(pmax-pmin),'Color',gray,'LineWidth',lw)
plot(Sim.tvec,(M.pbar-pmin)/(pmax-pmin),'linewidth',2,'color',col(2,:))
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel({'Probability   ';'of spiking    '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike train
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
BarVar=M.nbar+M.nvar;
BarVar(BarVar>1)=1;
stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',sw)
stem(Sim.tvec,M.nbar,'Marker','none','Color',col(2,:),'LineWidth',sw)
axis([xs min(R.n) max(R.n)])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel({'Spike        '; 'Train         '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike history
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,P.omega*S.h);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-hmin)/(hmax-hmin),ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(P.omega*R.h-hmin)/(hmax-hmin),'Color',gray,'LineWidth',lw)%(R.h-hmin)/(hmax-hmin)
plot(Sim.tvec,(P.omega*M.hbar-hmin)/(hmax-hmin),'linewidth',2,'color',col(2,:))
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel({'Spike History ';'Term        '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% calcium
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,S.C);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-cmin)/(cmax-cmin),ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray,'LineWidth',lw)%(R.C-cmin)/(cmax-cmin)
plot(Sim.tvec,(M.Cbar-cmin)/(cmax-cmin),'linewidth',2,'color',col(2,:))
% set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs 0 1])
ylab=ylabel({'Calcium     '; 'Concentration '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})

% observations
% subplot(Nrows,Ncols,12), cla, hold on, box on
% plot(Oind*Sim.dt,O,'.k','LineWidth',1,'markersize',7)
% set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
% axis([xs min(O(2:end-6)) max(O(2:end-6))])
% set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
text(-0.47, 5.9,'True','fontsize',tfs)
text(+0.18, 5.9,'Inferred','fontsize',tfs)

% make eps
fig=figure(1);
wh=[7 7];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\schem3