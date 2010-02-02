% this script does makes a figure showing how the smc-em approach
% can use stimulus information to further refine the stimulus
%
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nStimulus Simulation Fig\n')

%% 1) set simulation metadata

Sim.T       = 230;                                  % # of time steps
Sim.dt      = 1/100;                                % time step size
Sim.freq    = 5;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 100;                                  % # of particles
Sim.M       = 1;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = 1*randn(1,Sim.T);                     % stimulus

Sim.Mstep   = false;                                % do M-step
Sim.C_params = true;                                % whether to estimate calcium parameters {tau,A,C_0,sig}
% Sim.Rparams = true;                               % whether to estimate rate governing parameters {b,k}
% Sim.SpikHis = false;                              % whether to estimate spike history parameters {h}
% Sim.Oparams = false;                              % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 0;                                    % max # of EM iterartions

%% 2) initialize parameters

% initialize barrier and wiener filter parameters
P.rate  = 50/(Sim.T*Sim.dt);                        % expected spike rate
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 0.5;                                      % calcium decay time constant (sec)
P.lam   = Sim.T/(P.rate*P.A)*Sim.dt;                % expected jump size ber time bin
P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-P.rate/Sim.T)/Sim.dt);     % linear filter
P.tau_c     = P.tau;
P.A         = 15;
P.C_0       = 5;                                    % baseline [Ca++]
P.C_init    = P.C_0;                                % initial [Ca++]
P.sigma_c   = P.sig;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 100;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = 1e-5;                                 % scaled variance
P.zeta      = 4*P.gamma;                            % constant variance
P.a         = Sim.dt/P.tau_c;
            
if Sim.M==1                                         % if there are spike history terms
    P.omega = -1;                                   % weight
    P.tau_h = 0.02;                                 % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

%% 3) simulate data

kx        = P.k'*Sim.x;                                 %external input to neuron
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
U_sampl   = rand(1,Sim.T);                              %generate random number to use for sampling
if Sim.M==1                                              %if spike history terms, recursively
    p         = zeros(1,Sim.T);                       %extize p_t because it must be updated iteratively
    n         = zeros(1,Sim.T);                       %extize n_t because it must be updated iteratively
    h         = zeros(1,Sim.T);                   %extize spike history because it must be updated iteratively
    epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(1,Sim.T); %generate noise on spike history
    for t=2:Sim.T                                       %update states
        h(:,t)= (1-Sim.dt./P.tau_h).*h(:,t-1)+n(t-1) + epsilon_h(:,t);%update h terms
        y_t=kx(t)+P.omega'*h(:,t);                    %generate operand for rate function
        p(t)=1-exp(-exp(y_t)*Sim.dt);                 %generate rate
        n(t)  = U_sampl(t)<p(t);                    %sample from bernoulli with prob p_t
    end %time loop
else
    n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
end
C=zeros(1,Sim.T);
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
for t=2:Sim.T                                           %recursively update calcium
    C(t)  = (1-Sim.dt/P.tau_c)*C(t-1) + P.a*P.C_0 + P.A*n(t) + epsilon_c(t);
end
S=Hill_v1(P,C);
F=(P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(1,Sim.T))';
F(F<=0)=eps;
subplot(211), plot(F+1); hold on, stem(n),
if Sim.M==1, plot(P.omega*h); end
subplot(212),
if Sim.M==1, plot(p), hold on, end
plot(z1(Sim.x),'k');
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;          % for plotting purposes in ParticleFiltD

%% 4) infer spikes and estimate parameters

Algs=1:2;                                       % which algorithms within DataComp to use
for m=Algs
    Tim=Sim;
    if m==1, Tim.M=0; Tim.Alg=7; Tim.x=ones(size(Sim.x));
    elseif m==2, Tim.M=1; Tim.Alg=7; end
    I{m}    = DataComp13(F,P,Tim);
end
save StimSim

%% 5) plot results

load('StimSim.mat')
% Algs=[7 9];
fig     = figure(2); clf,
nrows   = numel(Algs)+6;
gray    = [.5 .5 .5];            % define gray color
col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter   = 'tex';                    % interpreter for axis labels
xlims   = [45 Sim.T-2*Sim.freq];   % xmin and xmax for current plot
fs      = 12;                       % font size
ms      = 4;                       % marker size for real spike
sw      = 2;                      % spike width
lw      = 2;                        % line width
xticks  = 0:20:200;               % XTick positions
% I{2}.name=[{'Wiener Filter'}];
I{1}.name=[{'Superresolution'}; {'PFS Spike Inference'}];
I{2}.name=[{'GLM PFS'}; {'Spike Inference'}];
spt = find(Sim.n==1);
tvec_o=xlims(1):Sim.freq:xlims(2);
i=0;

% plot stimulus
i=i+1; subplot(nrows,1,i), hold on
plot(Sim.x,'k','LineWidth',2);
ylab=ylabel([{'External'}; {'Stimulus'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims min(Sim.x) max(Sim.x)])

% plot spike history
i=i+1; subplot(nrows,1,i), hold on
plot(P.omega*h,'Color',gray,'LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Spike History'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims min(P.omega*h) 0.1])

% plot prob spiking
i=i+1; subplot(nrows,1,i), hold on
plot(p,'Color',gray,'LineWidth',2);
ylab=ylabel([{'Probability'}; {'of Spiking'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1])

% plot spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,Sim.n(spt),'Marker','none','LineWidth',sw,'Color',gray,'MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'Simulated'}; {'Spike Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1])

% plot calcium
i=i+1; subplot(nrows,1,i), hold on
plot(C/P.k_d,'Color',gray,'LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Calcium'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1*max(C/P.k_d)])

% plot fluorescence data
i=i+1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(F(tvec_o)),'.-k','LineWidth',2,'MarkerSize',15);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Simulated'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,
    BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
    spts=find(BarVar>1e-3);
    stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
    spts=find(I{m}.M.nbar>1e-3);
    stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
    stem(spt,Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
    axis([xlims 0 1])
    hold off,
    ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    box off
end

% label last subplot
set(gca,'XTick',tvec_o,'XTickLabel',[0;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''};.5;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''};1;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''};1.5;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''}])
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 7];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','StimSim')

%% 6) plot b&w results

load('StimSim.mat')
% Algs=[7 9];
fig     = figure(2); clf,
nrows   = numel(Algs)+6;
gray    = [.5 .5 .5];            % define gray color
col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter   = 'tex';                    % interpreter for axis labels
xlims   = [45 Sim.T-2*Sim.freq];   % xmin and xmax for current plot
fs      = 12;                       % font size
ms      = 4;                       % marker size for real spike
sw      = 2;                      % spike width
lw      = 2;                        % line width
xticks  = 0:20:200;               % XTick positions
% I{2}.name=[{'Wiener Filter'}];
I{1}.name=[{'Superresolution'}; {'PFS Spike Inference'}];
I{2}.name=[{'GLM PFS'}; {'Spike Inference'}];
spt = find(Sim.n==1);
tvec_o=xlims(1):Sim.freq:xlims(2);
i=0;

% plot stimulus
i=i+1; subplot(nrows,1,i), hold on
plot(z1(Sim.x),'k','LineWidth',2);
ylab=ylabel([{'External'}; {'Stimulus'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot spike history
i=i+1; subplot(nrows,1,i), hold on
plot(z1(P.omega*h)-1,'Color','k','LineWidth',2,'LineStyle','-');
ylab=ylabel([{'Simulated'}; {'Spike History'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',-1,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims -1.1 0])

% plot prob spiking
i=i+1; subplot(nrows,1,i), hold on
plot(p,'Color','k','LineWidth',2);
ylab=ylabel([{'Probability'}; {'of Spiking'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,Sim.n(spt),'Marker','none','LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'Simulated'}; {'Spike Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1])

% plot calcium
i=i+1; subplot(nrows,1,i), hold on
plot(z1(C/P.k_d),'Color','k','LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Calcium'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot fluorescence data
i=i+1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(F(tvec_o)),'.-k','LineWidth',2,'MarkerSize',15);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Simulated'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',1,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot inferred spike trains
for m=Algs
    i=i+1;
    subplot(nrows,1,i), hold on,
    BarVar=I{m}.M.nbar+I{m}.M.nvar; BarVar(BarVar>1)=1;
    spts=find(BarVar>1e-3);
    stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',gray);
    spts=find(I{m}.M.nbar>1e-3);
    stem(spts,I{m}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color','k')
    stem(spt,1.1*Sim.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
    axis([xlims 0 1.1])
    hold off,
    ylab=ylabel(I{m}.name,'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',1,'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    box off
end

% label last subplot
set(gca,'XTick',tvec_o,'XTickLabel',[0;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''};.5;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''};1;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''};1.5;{''};{''};{''};{''};{''};{''};{''};{''};{''};{''}])
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 7];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','StimSim_bw')