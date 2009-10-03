% this script does makes a figure showing how the smc-em approach
% outperforms other approaches in saturating data:
% 
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nSaturation Simulation Fig\n')

%% 1) set simulation metadata

Sim.T       = 405;                                  % # of time steps
Sim.dt      = 0.025;                                % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 200;                                  % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus
Sim.rate    = 20;                                   % expected spike rate

Sim.Mstep   = false;                                % do M-step
Sim.C_params = true;                                % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = true;                                % whether to estimate rate governing parameters {b,k}
Sim.h_params = false;                               % whether to estimate spike history parameters {h}
Sim.F_params = false;                               % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 0;                                    % max # of EM iterartions

%% 2) initialize parameters

% initialize barrier and wiener filter parameters
P.rate  = 10/(Sim.T*Sim.dt);
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 2;                                        % calcium decay time constant (sec)
P.lam   = Sim.T/(Sim.rate*P.A)*Sim.dt;              % expected jump size ber time bin
P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-Sim.rate/Sim.T)/Sim.dt);   % linear filter
P.tau_c     = P.tau;
P.A         = 50;
P.C_0       = 20;                                   % baseline [Ca++]
P.C_init    = P.C_0;                                % initial [Ca++]
P.sigma_c   = P.sig;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 1e-4;                                 % scaled variance
P.zeta      = 4*P.gamma;                            % constant variance
P.a         = Sim.dt/P.tau_c;

%% 3) simulate data

n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
C=zeros(Sim.T,1);
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
for t=2:Sim.T                                           %recursively update calcium
    C(t)  = (1-P.a)*C(t-1) + P.a*P.C_0 + P.A*n(t) + epsilon_c(t);   
end
S=Hill_v1(P,C);
F=P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(Sim.T,1);
F(F<=0)=eps;
figure(1), plot(F+1); hold on, stem(n)
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;                 % for plotting purposes in ParticleFiltD

%% 4) infer spikes and estimate parameters

Algs=[2 7];                                       % which algorithms within DataComp to use
Sim.freq = 1; 
for m=Algs
    Sim.Alg = m;
    I{m}    = DataComp13(F,P,Sim);
end
save SaturSim


%% 5) plot results

load('SaturSim.mat')

fig=figure(2); clf,
nrows = 6;
gray  = [.5 .5 .5];                 % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'tex';                      % interpreter for axis labels
xlims = [0+2 Sim.T-4*Sim.freq];     % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                               % marker size for real spike
sw=1.5;                             % spike width
lw=2;                               % line width
% make xticks
Nsec = floor(Sim.T*Sim.dt);
secs = zeros(1,Nsec-1);
for i=1:Nsec
    secs(i) = find(Sim.tvec>=i,1);
end
I{2}.name=[{'Wiener Filter'}];
I{7}.name=[{'NOOPSI Filter'}; {'Spike Inference'}];
spt=find(Sim.n==1);
i=0;

col   = [1 0 0; 0.2 0.2 1];             % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;         % define colors for std
% I{7}.name=[{'Linear Observation'}; {'Particle Filter'}];
I{7}.name=[{'Nonlinear Observation'}; {'PFS Spike Inference'}];
% I{7}.name=[{'Superresolution'}; {'Particle Filter'}];
% I{7}.name=[{'GLM'}; {'Particle Filter'}];


% plot spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,Sim.n(spt),'Marker','none','LineWidth',sw,'Color',gray,'MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'Simulated'}; {'Spike Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[],'XTickLabel',[])
axis([xlims 0 1])

% plot calcium
i=i+1; subplot(nrows,1,i), hold on
plot(C/P.k_d,'Color',gray,'LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Calcium'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTickLabel',[])
axis([xlims 0 1.1*max(C/P.k_d)])

% plot fluorescence data
i=i+1; subplot(nrows,1,i), hold on
plot(z1(F),'k','LineWidth',2);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Simulated'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[],'XTickLabel',[])
axis([xlims 0 1])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if m==8,
        plot(D{j}.F+1,'.','Color',gray)
        plot(I{m}.P.alpha*Hill_v1(I{m}.P,I{m}.C)+1,'Color',col(2,:),'linewidth',2),
        stem(SubSampleSpikes(D,freq),'Marker','o','LineWidth',sw,'Color',gray)
        stem(I{m}.M.nbar,'Marker','none','LineWidth',sw,'Color',col(2,:))
        axis([xlims 0 1])
    elseif any(m==[7 9])
        BarVar=I{m}.M.nbar+(I{m}.M.nvar); BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        stem(spt,Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
        axis([xlims 0 1])
    else
        n_est   = I{m}.n; n_est   = n_est/max(n_est);   %normalize estimate
        neg = find(n_est<=0);
        stem(neg,n_est(neg),'Marker','none','LineWidth',sw,'Color',col(1,:))
        pos = find(n_est>0);
        stem(pos,n_est(pos),'Marker','none','LineWidth',sw,'Color',col(2,:))
        stem(spt,Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')

        axis([xlims min(n_est) max(n_est)])
    end

    hold off,
    ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',secs,'XTickLabel',[])
end

% plot calcium
if any(Algs==7) && ~any(Algs==9)
    subplot(nrows,1,i+1), hold on
%     plot(C/P.k_d,'k','LineWidth',2);
    C_est = I{m}.M.Cbar/I{m}.P.k_d;
    hfill=fill([1:Sim.T Sim.T:-1:1],[I{m}.M.Cptiles(.05*Sim.N,:) I{m}.M.Cptiles(.95*Sim.N,Sim.T:-1:1)]/I{m}.P.k_d,ccol(2,:));
    set(hfill,'edgecolor',ccol(2,:))
    plot(C_est,'Color',col(2,:),'LineWidth',2)
    set(gca,'YTick',1,'YTickLabel',[])
    set(gca,'XTick',secs,'XTickLabel',[])
    ylab=ylabel([{'Nonlinear Observation'}; {'PFS [Ca^2^+] Inference'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis([xlims 0 1.1*max(C_est)])
end

% plot stimulus
if any(Algs==9)
    subplot(nrows,1,i+1), hold on
    plot(x*10^-3,'Color','k','LineWidth',lw)
    ylab=ylabel([{'Stimulus'};{'(\muA)'}],'Interpreter',inter,'FontSize',fs,'FontName','Arial');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'XTick',secs,'XTickLabel',[])
    axis([xlims min(x)*10^-3 max(x)*10^-3])
end
set(gca,'XTick',secs,'XTickLabel',1:Nsec,'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 6];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SaturSim')

%% 6) plot b&w results

load('SaturSim.mat')

fig=figure(2); clf,
nrows = 6;
gray  = [.5 .5 .5];                 % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'tex';                      % interpreter for axis labels
xlims = [0+2 Sim.T-4*Sim.freq];     % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                               % marker size for real spike
sw=1.5;                             % spike width
lw=2;                               % line width
% make xticks
Nsec = floor(Sim.T*Sim.dt);
secs = zeros(1,Nsec-1);
for i=1:Nsec
    secs(i) = find(Sim.tvec>=i,1);
end
I{2}.name=[{'Wiener Filter'}];
I{7}.name=[{'NOOPSI Filter'}; {'Spike Inference'}];
spt=find(Sim.n==1);
i=0;

col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
% I{7}.name=[{'Linear Observation'}; {'Particle Filter'}];
I{7}.name=[{'Nonlinear Observation'}; {'PFS Spike Inference'}];
% I{7}.name=[{'Superresolution'}; {'Particle Filter'}];
% I{7}.name=[{'GLM'}; {'Particle Filter'}];


% plot spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,Sim.n(spt),'Marker','none','LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'Simulated'}; {'Spike Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[],'XTickLabel',[])
axis([xlims 0 1])

% plot calcium
i=i+1; subplot(nrows,1,i), hold on
plot(C/P.k_d,'Color','k','LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Calcium'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTickLabel',[])
axis([xlims 0 1.1*max(C/P.k_d)])

% plot fluorescence data
i=i+1; subplot(nrows,1,i), hold on
plot(z1(F),'k','LineWidth',2);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Simulated'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[],'XTickLabel',[])
axis([xlims 0 1])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if m==8,
        plot(D{j}.F+1,'.','Color',gray)
        plot(I{m}.P.alpha*Hill_v1(I{m}.P,I{m}.C)+1,'Color',col(2,:),'linewidth',2),
        stem(SubSampleSpikes(D,freq),'Marker','o','LineWidth',sw,'Color',gray)
        stem(I{m}.M.nbar,'Marker','none','LineWidth',sw,'Color',col(2,:))
        axis([xlims 0 1])
    elseif any(m==[7 9])
        BarVar=I{m}.M.nbar+(I{m}.M.nvar); BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',gray);
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color','k')
        stem(spt,1.1*Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
        axis([xlims 0 1.1])
    else
        n_est   = I{m}.n; n_est   = n_est/max(n_est);   %normalize estimate
        neg = find(n_est<=0);
        stem(neg,n_est(neg),'Marker','none','LineWidth',sw,'Color',gray)
        pos = find(n_est>0);
        stem(pos,n_est(pos),'Marker','none','LineWidth',sw,'Color','k')
        stem(spt,1.1*Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)

        axis([xlims min(n_est(xlims(1):xlims(2))) 1.1])
    end

    hold off,
    ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',secs,'XTickLabel',[])
end

% plot calcium
if any(Algs==7) && ~any(Algs==9)
    subplot(nrows,1,i+1), hold on
%     plot(C/P.k_d,'k','LineWidth',2);
    C_est = I{m}.M.Cbar/I{m}.P.k_d;
    hfill=fill([1:Sim.T Sim.T:-1:1],[I{m}.M.Cptiles(.05*Sim.N,:) I{m}.M.Cptiles(.95*Sim.N,Sim.T:-1:1)]/I{m}.P.k_d,gray);
    set(hfill,'edgecolor',gray)
    plot(C_est,'Color','k','LineWidth',2)
%     plot(C/P.k_d,'Color','w','LineWidth',.5,'LineStyle','-.');
    set(gca,'YTick',1,'YTickLabel',[])
    set(gca,'XTick',secs,'XTickLabel',[])
    ylab=ylabel([{'Nonlinear Observation'}; {'PFS [Ca^2^+] Inference'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis([xlims 0 1.1*max(C_est)])
end

% plot stimulus
if any(Algs==9)
    subplot(nrows,1,i+1), hold on
    plot(x*10^-3,'Color','k','LineWidth',lw)
    ylab=ylabel([{'Stimulus'};{'(\muA)'}],'Interpreter',inter,'FontSize',fs,'FontName','Arial');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'XTick',secs,'XTickLabel',[])
    axis([xlims min(x)*10^-3 max(x)*10^-3])
end
set(gca,'XTick',secs,'XTickLabel',1:Nsec,'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 6];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SaturSim_bw')