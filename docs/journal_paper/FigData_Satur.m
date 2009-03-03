% this script does makes a figure showing how the smc-em approach
% outperforms linear algorithms with an example data set by doing:
%
% 1) load some data and get some parameters using real data
% 2) set simulation metadata (eg, dt, T, # particles, etc.)
% 3) initialize parameters
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nSaturation Fig\n')

%% 1) get data

% load '/Users/joshyv/Research/data/rafa_data/last_2_weeks/080218a/D080218a.mat';
load '/Users/joshyv/Research/data/rafa/brendon/last_2_weeks/080218a/D080218a.mat';
j=24;                                               % experiment number

% arrange stuff for plotting purposes
freq        = 1;
D{j}.spt    = Get_spt(D{j}.V);                      % get spike times
D{j}.n      = SubSampleSpikes(D{j},freq);              % make spike train
T           = min(length(D{j}.F),length(D{j}.n));   % # of time steps
Sim.n       = D{j}.n; %[D{j}.n; zeros(T-numel(D{j}.n),1)]; %Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD
F           = NaN*ones(T,1);
tvec_o      = freq:freq:T;
F(tvec_o)   = D{j}.F(1:T);
figure(1), clf, plot(F,'.'), hold on, stem(Sim.n,'k')

% [x C F] = GetDataParams(D{j}); %x0 = [A a C_0 k_d F_max];

%% 2) set simulation metadata

Sim.T_o     = T;                                    % # of time steps
Sim.dt      = D{j}.dt_o;                            % time step size
Sim.freq    = freq;                                    % # of time steps between observations
Sim.T       = Sim.T_o*Sim.freq;                     % # of observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 100;                                  % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus

Sim.Mstep   = true;                                 % do M-step
Sim.C_params = true;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params = false;                                % whether to estimate spike history parameters {h}
Sim.F_params = false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 50;                                    % max # of EM iterartions

%% 3) initialize parameters

% initialize barrier and wiener filter parameters
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 2;                                        % calcium decay time constant (sec)
P.lam   = Sim.T/(sum(Sim.n(1:Sim.T))*P.A)*Sim.dt;  % expected jump size ber time bin
P.sig   = 0.1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-sum(Sim.n(1:Sim.T))/Sim.T)/Sim.dt);% linear filter
P.tau_c     = 1;%Sim.dt/x(2);
P.A         = 20;
P.C_0       = 25;                                   % baseline [Ca++]
P.C_init    = 0.02;                                 % initial [Ca++]
P.sigma_c   = 10;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 2;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 1e-5;                                 % scaled variance
P.zeta      = 5e-5;                                 % constant variance

%% 4) infer spikes and estimate parameters

Algs=[2 7];                                       % which algorithms within DataComp to use
Sim.freq = 1;
for m=Algs
    Sim.Alg = m;
    if m==6, Sim.Alg=7; Sim.pf=0; end
    I{m}    = DataComp14(D{j}.F(1:Sim.T),P,Sim);
end
% save SaturData

%% 5) plot results
load('SaturData.mat')

fig=figure(2); clf,
if any(Algs==7) || any(Algs==9), nrows = numel(Algs)+2;
else nrows = numel(Algs)+1; end
gray  = [.75 .75 .75];              % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'tex';                     % interpreter for axis labels
xlims = [70 Sim.T-2*Sim.freq];      % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                              % marker size for real spike
sw=1.5;                             % spike width
lw=2;                               % line width
I{2}.name=[{'Wiener Filter'}];
I{4}.name=[{'Optimal'}; {'Non-negative Filter'}];
I{7}.name=[{'Nonlinear Observation'}; {'PFS Spike Inference'}];
% I{9}.name=[{'Superresolution'}; {'Particle Filter'}];
spt=find(D{j}.n==1); Sim.n=zeros(size(Sim.n)); Sim.n(spt)=1;
% make xticks
Nsec = floor(D{j}.T_o*D{j}.dt_o);
secs = zeros(1,Nsec-2);
for i=1:Nsec-2
    secs(i) = find(D{j}.FrameStartTime>i,1);
end

col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std

% plot real data
i=1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(D{j}.F(tvec_o)),'-k','LineWidth',2,'MarkerSize',ms);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'in vitro'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 1.1*max(D{j}.F(1:xlims(2)))])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if any(m==[6 7 9])
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        stem(spt,Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
        axis([xlims 0 1])
    else
        n_est   = I{m}.n; n_est   = n_est/max(n_est);   %normalize estimate

        neg = find(n_est<=0);
        stem(neg,n_est(neg),'Marker','none','LineWidth',sw,'Color',col(1,:))
        pos = find(n_est>0);
        stem(pos,n_est(pos),'Marker','none','LineWidth',sw,'Color',col(2,:))

        stem(spt,Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)

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
    C = I{m}.M.Cbar/I{m}.P.k_d;
    hfill=fill([1:Sim.T Sim.T:-1:1],[I{m}.M.Cptiles(1,:) I{m}.M.Cptiles(end,Sim.T:-1:1)]/I{m}.P.k_d,ccol(2,:));
    set(hfill,'edgecolor',ccol(2,:))
    plot(C(1:end-1),'Color',col(2,:),'LineWidth',2)
    set(gca,'YTick',1,'YTickLabel',[])
    ylab=ylabel([{'Nonlinear Observation'}; {'PFS [Ca^{2+}] Inference'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis([xlims min(C) 1.1*max(C(1:xlims(2)))])
end

% plot stimulus
if any(Algs==9)
    subplot(nrows,1,i+1), hold on
    plot(Sim.x(2,:),'Color','k','LineWidth',lw)
    ylab=ylabel([{'Stimulus'};{'(mA)'}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[-.5 .5],'YTickLabel',round([Sim.minx/2 Sim.maxx/2]*1e-3))
    axis([xlims min(Sim.x(2,:)) max(Sim.x(2,:))])
end

set(gca,'XTick',secs,'XTickLabel',round((secs-xlims(1))*Sim.dt),'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 4];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SaturData')

%% 5) plot b&w results
load('SaturData.mat')

fig=figure(2); clf,
if any(Algs==7) || any(Algs==9), nrows = numel(Algs)+2;
else nrows = numel(Algs)+1; end
gray  = [.75 .75 .75];              % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'tex';                     % interpreter for axis labels
xlims = [70 Sim.T-2*Sim.freq];      % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                              % marker size for real spike
sw=1.5;                             % spike width
lw=2;                               % line width
I{2}.name=[{'Wiener Filter'}];
I{4}.name=[{'Optimal'}; {'Non-negative Filter'}];
I{7}.name=[{'Nonlinear Observation'}; {'PFS Spike Inference'}];
% I{9}.name=[{'Superresolution'}; {'Particle Filter'}];
spt=find(D{j}.n==1); Sim.n=zeros(size(Sim.n)); Sim.n(spt)=1;
% make xticks
Nsec = floor(D{j}.T_o*D{j}.dt_o);
secs = zeros(1,Nsec-2);
for i=1:Nsec-2
    secs(i) = find(D{j}.FrameStartTime>i,1);
end

col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std

% plot real data
i=1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(D{j}.F(tvec_o)),'-k','LineWidth',2,'MarkerSize',ms);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'in vitro'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 1.1*max(D{j}.F(1:xlims(2)))])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if any(m==[6 7 9])
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
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
    C = I{m}.M.Cbar/I{m}.P.k_d;
    hfill=fill([1:Sim.T Sim.T:-1:1],[I{m}.M.Cptiles(1,:) I{m}.M.Cptiles(end,Sim.T:-1:1)]/I{m}.P.k_d,ccol(2,:));
    set(hfill,'edgecolor',gray)
    plot(C(1:end-1),'Color','k','LineWidth',2)
    set(gca,'YTick',1,'YTickLabel',[])
    ylab=ylabel([{'Nonlinear Observation'}; {'PFS [Ca^{2+}] Inference'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis([xlims min(C) 1.1])
end

% plot stimulus
if any(Algs==9)
    subplot(nrows,1,i+1), hold on
    plot(Sim.x(2,:),'Color','k','LineWidth',lw)
    ylab=ylabel([{'Stimulus'};{'(mA)'}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[-.5 .5],'YTickLabel',round([Sim.minx/2 Sim.maxx/2]*1e-3))
    axis([xlims min(Sim.x(2,:)) max(Sim.x(2,:))])
end

set(gca,'XTick',secs,'XTickLabel',round((secs-xlims(1))*Sim.dt),'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 4];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','SaturData_bw')