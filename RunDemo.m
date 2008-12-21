% this script does makes a figure showing how the smc-em approach
% can use stimulus information to further refine the stimulus
%
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nDemo\n')

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.T       = 230;                                  % # of time steps
Sim.dt      = 1/100;                                % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 100;                                  % # of particles
Sim.M       = 1;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = 1*randn(1,Sim.T);                     % stimulus

Sim.Mstep   = 1;                                    % do M-step
Sim.C_params= 1;                                    % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= 0;                                    % whether to estimate rate governing parameters {b,k}
Sim.h_params= 0;                                    % whether to estimate spike history parameters {h}
Sim.F_params= 0;                                    % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 10;                                    % max # of EM iterartions

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kx        = P.k'*Sim.x;                             %external input to neuron
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);  %generate noise on calcium
U_sampl   = rand(1,Sim.T);                          %generate random number to use for sampling
if Sim.M==1                                         %if spike history terms, recursively
    p         = zeros(1,Sim.T);                     %extize p_t because it must be updated iteratively
    n         = zeros(1,Sim.T);                     %extize n_t because it must be updated iteratively
    h         = zeros(1,Sim.T);                     %extize spike history because it must be updated iteratively
    epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(1,Sim.T); %generate noise on spike history
    for t=2:Sim.T                                   %update states
        h(:,t)= (1-Sim.dt./P.tau_h).*h(:,t-1)+n(t-1) + epsilon_h(:,t);%update h terms
        y_t=kx(t)+P.omega'*h(:,t);                  %generate operand for rate function
        p(t)=1-exp(-exp(y_t)*Sim.dt);               %generate rate
        n(t)  = U_sampl(t)<p(t);                    %sample from bernoulli with prob p_t
    end %time loop
else
    n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
end
C=zeros(1,Sim.T);
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);  %generate noise on calcium
for t=2:Sim.T                                       %recursively update calcium
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
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;             % for plotting purposes in ParticleFiltD

%% 4) infer spikes and estimate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

P.A = 2*P.A;
[I.M I.P]   = GOOPSI_main_v1_0(F,P,Sim);
I.n         = I.M.nbar;
save StimSim

%% 5) plot results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('StimSim.mat')
fig     = figure(2); clf,
nrows   = 7;
gray    = [.75 .75 .75];            % define gray color
col   = [1 0 0; 0.2 0.2 1];         % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter   = 'tex';                    % interpreter for axis labels
xlims   = [45 Sim.T-2*Sim.freq];    % xmin and xmax for current plot
fs      = 12;                       % font size
ms      = 20;                       % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
xticks  = 0:20:Sim.T;               % XTick positions
spt = find(Sim.n==1);               % find spike times
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
plot(P.omega*h,'k','LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Spike History'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims min(P.omega*h) 0.1])

% plot prob spiking
i=i+1; subplot(nrows,1,i), hold on
plot(p,'k','LineWidth',2);
ylab=ylabel([{'Probability'}; {'of Spiking'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1])

% plot spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'Simulated'}; {'Spike Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1])

% plot calcium
i=i+1; subplot(nrows,1,i), hold on
plot(C/P.k_d,'k','LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Calcium'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1*max(C/P.k_d)])

% plot fluorescence data
i=i+1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(F(tvec_o)),'.-k','LineWidth',2,'MarkerSize',ms*.75);
% stem(Sim.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Simulated'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot inferred spike trains
i=i+1; subplot(nrows,1,i), hold on,
stem(spt,Sim.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
BarVar=I.M.nbar+I.M.nvar; BarVar(BarVar>1)=1;
spts=find(BarVar>1e-3);
stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
spts=find(I.M.nbar>1e-3);
stem(spts,I.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
axis([xlims 0 1])
hold off,
ylab=ylabel([{'Particle'}; {'Filter'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
box off

% label last subplot
set(gca,'XTick',xticks,'XTickLabel',(xticks-xlims(1))*Sim.dt)
xlabel('Time (sec)','FontSize',fs)

% print fig
wh=[7 7];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','StimSim')