% this script does makes a figure showing how the smc-em approach
% outperforms can incorporate prior information to refine the spike train
% estimates by doing the following:
%
% 1) load some data
% 2) set simulation metadata (eg, dt, T, # particles, etc.)
% 3) initialize parameters
% 4) infers spikes using a variety of approaches
% 5) plots results
clear, clc, close all, fprintf('\nStimulus Fig\n')

%% 1) get data

load '/Users/joshyv/Research/data/rafa/brendon/last_2_weeks/080218a/D080218a.mat';
j=24;                                               % experiment number

D{j}.spt    = Get_spt(D{j}.V);                      % get spike times
D{j}.n      = SubSampleSpikes(D{j},1);              % make spike train

%% 2) set simulation metadata

Sim.T       = min(length(D{j}.n),length(D{j}.F));   % # of time steps
Sim.dt      = D{j}.dt_o;                            % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 200;                                  % # of particles
Sim.M       = 1;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.n       = D{j}.n; Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD

Sim.Mstep   = true;                                 % do M-step
Sim.C_params= true;                                % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params= true;                                % whether to estimate spike history parameters {h}
Sim.F_params= false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 50;                                   % max # of EM iterartions

%% 3) initialize parameters

% initialize barrier and wiener filter parameters
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 2;                                        % calcium decay time constant (sec)
P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.tau_c     = 1.2;
P.A         = 16.5;
P.C_0       = 16;                                   % baseline [Ca++]
P.C_init    = 0.02;                                 % initial [Ca++]
P.sigma_c   = 10;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 2;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 2e-5;                                 % scaled variance
P.zeta      = 5*P.gamma;                                 % constant variance

if Sim.M==1                                         % if there are spike history terms
    P.omega = -0.5;                                   % weight
    P.tau_h = 0.015;                                    % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

%% 4) infer spikes and estimate parameters

Algs=[2 7 9];                                           % which algorithms within DataComp to use
Sim.freq    = 4;
tvec_o      = Sim.freq:Sim.freq:Sim.T;              % observation time steps
for m=Algs
    Sim.Alg = m;
    if m==2,
        F       = D{j}.F(tvec_o);
        Tim     = Sim;
        Tim.dt  = Sim.dt*Sim.freq;
        Tim.T   = Sim.T/Sim.freq;
        P.lam   = Sim.T/(sum(D{j}.n(1:Sim.T))*P.A)*Sim.dt;  % expected jump size ber time bin
    elseif m==6 %do prior sampler
        F           = D{j}.F(1:Sim.T);
        Sim.StimDim = 1;                                    % # of stimulus dimensions
        Sim.x       = ones(Sim.StimDim,Sim.T);
        Sim.Alg     = 7;
        Sim.pf      = 0;
        Tim         = Sim;
        P.k         = [log(-log(1-sum(D{j}.n)/Sim.T)/Sim.dt)];% linear filter
    elseif m==7
        F           = D{j}.F(1:Sim.T);
        Sim.StimDim = 1;                                    % # of stimulus dimensions
        Sim.x       = ones(Sim.StimDim,Sim.T);
        Tim         = Sim;
        P.k         = [log(-log(1-sum(D{j}.n)/Sim.T)/Sim.dt)];% linear filter
    elseif m==9
        F           = D{j}.F(1:Sim.T);
        Sim.Alg     = 7;
        Sim.StimDim = 2;                                    % # of stimulus dimensions
        x           = SubSampleStim(D{j},1)';               % incorporate stimulus
        Sim.x       = ones(Sim.StimDim,numel(x));
        Sim.maxx    = max(x);
        Sim.minx    = min(x);
        Sim.x(2,:)  = x/Sim.maxx;                     % normalize stimuls
        Tim         = Sim;
        P.k         = [log(-log(1-sum(D{j}.n)/Sim.T)/Sim.dt); 1];% linear filter
    end
    I{m} = DataComp13(F,P,Tim);
end
save StimData2

%% 5) plot results

load('StimData2.mat')

fig=figure(2); clf,
if any(Algs==7) || any(Algs==9), nrows = numel(Algs)+3;
else nrows = numel(Algs)+1; end
gray  = [.5 .5 .5];              % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'tex';                     % interpreter for axis labels
xlims = [80 Sim.T-3];               % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                              % marker size for real spike
sw=2.5;                             % spike width
sw2=8;                              % spike width for wiener
lw=2;                               % line width
I{7}.name=[{'Superresolution'}; {'PFS Spike Inference'}];
I{9}.name=[{'GLM PFS'};  {'Spike Inference'}];
Sim.n=D{j}.n;
spt=find(Sim.n==1);
% make xticks
Nsec = floor(D{j}.T_o*D{j}.dt_o);
secs = zeros(1,Nsec-1);
for i=1:Nsec
    secs(i) = find(D{j}.FrameStartTime>i,1);
end
i=0;

col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std


% plot stimulus
if any(Algs==9)
    i=i+1; subplot(nrows,1,i), hold on
    plot(Sim.x(2,:),'Color','k','LineWidth',lw)
    ylab=ylabel([{'Stimulus'}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
%     ylab=ylabel([{'Stimulus'};{['(\mu' 'A)']}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[-.5 .5],'YTickLabel',[]);%[{'-36'}; {'+35'}]); %round([Sim.minx/2 Sim.maxx/2]*1e-3))
    set(gca,'XTick',secs,'XTickLabel',[])
    axis([xlims min(Sim.x(2,:)) max(Sim.x(2,:))])
end

% plot real spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,Sim.n(spt),'Marker','none','LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'in vitro'}; {'Spikes'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 1])

% plot real fluorescence
i=i+1; subplot(nrows,1,i), hold on
plot(tvec_o,(D{j}.F(tvec_o)),'-','Color',gray,'LineWidth',2,'MarkerSize',15);
plot(tvec_o,(D{j}.F(tvec_o)),'.k','LineWidth',2,'MarkerSize',15);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'in vitro'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims min(D{j}.F(xlims(1):xlims(2))) max(D{j}.F(xlims(1):xlims(2)))])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if any(m==[6 7 9])
%         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        stem(spt,Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
        axis([xlims 0 1])
    else
%         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        plot(zeros(size(Sim.n)),'k')
        n_est   = I{m}.n; n_est   = n_est/max(n_est);   %normalize estimate
        
                
        neg = find(n_est<=0);
        stem(tvec_o(neg),n_est(neg),'Marker','none','LineWidth',sw2,'Color',col(1,:))        
        pos = find(n_est>0);
        stem(tvec_o(pos),n_est(pos),'Marker','none','LineWidth',sw2,'Color',col(2,:))

%         stem(Sim.freq:Sim.freq:Sim.T,n_est,'Marker','none','LineWidth',8,'Color',col(1,:));
        %         area(tvec_o,n_est,'FaceColor',col(1,:),'EdgeColor',col(1,:))
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
    C = I{m}.M.Cbar/I{m}.P.k_d;
    hfill=fill([1:Sim.T Sim.T:-1:1],[I{m}.M.Cptiles(1,:) I{m}.M.Cptiles(end,Sim.T:-1:1)]/I{m}.P.k_d,ccol(2,:));
    set(hfill,'edgecolor',ccol(2,:))
    plot(C,'Color',col(2,:),'LineWidth',2)
    set(gca,'YTick',1,'YTickLabel',[])
    ylab=ylabel([{'Inferred'};{'Calcium'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis([xlims min(C) max(C)+.1])
end

set(gca,'XTick',secs,'XTickLabel',round((secs-xlims(1))*Sim.dt),'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','StimData')

%% 6) plot b&w results

load('StimData2.mat')

fig=figure(2); clf,
if any(Algs==7) || any(Algs==9), nrows = numel(Algs)+3;
else nrows = numel(Algs)+1; end
gray  = [.5 .5 .5];              % define gray color
gray2  = [.7 .7 .7];              % define gray color
col     = [1 0 0; 0 .5 0];          % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   % define colors for std
inter = 'tex';                     % interpreter for axis labels
xlims = [80 Sim.T-3];               % xmin and xmax for current plot
fs=12;                              % font size
ms=4;                              % marker size for real spike
sw=2.5;                             % spike width
sw2=8;                              % spike width for wiener
lw=2;                               % line width
I{7}.name=[{'Superresolution'}; {'PFS Spike Inference'}];
I{9}.name=[{'GLM PFS'};  {'Spike Inference'}];
Sim.n=D{j}.n;
spt=find(Sim.n==1);
% make xticks
Nsec = floor(D{j}.T_o*D{j}.dt_o);
secs = zeros(1,Nsec-1);
for i=1:Nsec
    secs(i) = find(D{j}.FrameStartTime>i,1);
end
i=0;

col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std


% plot stimulus
if any(Algs==9)
    i=i+1; subplot(nrows,1,i), hold on
    plot(Sim.x(2,:),'Color','k','LineWidth',lw)
    ylab=ylabel([{'Stimulus'}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
%     ylab=ylabel([{'Stimulus'};{['(\mu' 'A)']}],'Interpreter',inter,'FontSize',fs,'FontName','Helvetica');
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[-.5 .5],'YTickLabel',[]);%[{'-36'}; {'+35'}]); %round([Sim.minx/2 Sim.maxx/2]*1e-3))
    set(gca,'XTick',secs,'XTickLabel',[])
    axis([xlims min(Sim.x(2,:)) max(Sim.x(2,:))])
end

% plot real spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,Sim.n(spt),'Marker','none','LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k','LineStyle','-')
ylab=ylabel([{'in vitro'}; {'Spikes'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims 0 1])

% plot real fluorescence
i=i+1; subplot(nrows,1,i), hold on
plot(tvec_o,(D{j}.F(tvec_o)),'-','Color',gray,'LineWidth',2,'MarkerSize',15);
plot(tvec_o,(D{j}.F(tvec_o)),'.k','LineWidth',2,'MarkerSize',15);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'in vitro'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',secs,'XTickLabel',[])
axis([xlims min(D{j}.F(xlims(1):xlims(2))) max(D{j}.F(xlims(1):xlims(2)))])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,

    if any(m==[6 7 9])
%         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',gray);
        spts=find(I{m}.M.nbar>1e-3);
        stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color','k')
        stem(spt,1.1*Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','Color','k','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
        axis([xlims 0 1.1])
    else
%         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        plot(zeros(size(Sim.n)),'k')
        n_est   = I{m}.n; n_est   = n_est/max(n_est);   %normalize estimate
                
        neg = find(n_est<=0);
        stem(tvec_o(neg),n_est(neg),'Marker','none','LineWidth',sw2,'Color',gray)        
        pos = find(n_est>0);
        stem(tvec_o(pos),n_est(pos),'Marker','none','LineWidth',sw2,'Color','k')

%         stem(Sim.freq:Sim.freq:Sim.T,n_est,'Marker','none','LineWidth',8,'Color',col(1,:));
        %         area(tvec_o,n_est,'FaceColor',col(1,:),'EdgeColor',col(1,:))
        stem(spt,1.1*Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
        axis([xlims min(n_est) 1.1])
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
    plot(C,'Color',col(2,:),'LineWidth',2)
    set(gca,'YTick',1,'YTickLabel',[])
    ylab=ylabel([{'Inferred'};{'Calcium'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis([xlims min(C) max(C)+.1])
end

set(gca,'XTick',secs,'XTickLabel',round((secs-xlims(1))*Sim.dt),'FontSize',fs)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 5];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','StimData_bw')