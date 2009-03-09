%This script fits a cell and prepares parameters for glm-SMC sampler
clear, clc, fprintf('\nfit cell ')

%CELL sample files
%make file contain fluorescence in var D{2}.F and the list of spike times in R.spt
% fname=['../spatial-filter/TS108_6_003.mat']; fprintf('%s\n',fname)
fname=['../spatial-filter/TS108_6_012.mat']; fprintf('%s\n',fname)
% fname=['../spatial-filter/FilteredData_D080218a24.mat']; fprintf('%s\n',fname)

clear D R
load(fname);


%Takashi's CELL
dt=0.02;                             % time step
F=0.99*D{2}.F+0.01;                  % fluorescence from data
n=zeros(size(F)); n(R.spt)=1;        % spike train from data

% %Rafa's CELL
% dt=0.026;
% n=D.n;
% F=0.99*D.F(1:length(n))+0.1;

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim.T       = length(F);             % # of time steps
Sim.dt      = dt;                    % time step size
Sim.freq    = 1;                     % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;          % # of actual seconds
Sim.T_o     = Sim.T;                 % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;% vector of times
Sim.N       = 25;                    % # of particles
Sim.M       = 0;                     % # of spike history dimensions
Sim.pf      = 1;                     % conditional sampler (not prior)?
Sim.StimDim = 1;                     % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);         % stimulus, uniform firing

Sim.Mstep   = 1;                     % do M-step
Sim.MaxIter = 15;                    % # of EM iterartions
Sim.C_params= 1;                     % estimate calcium model {tau,A,C_0,sig}
Sim.n_params= 1;                     % estimate rate model {b,k}
Sim.h_params= 1;                     % estimate spike history model {h}
Sim.F_params= 1;                     % estimate observation model {alpha,beta}
Sim.G_params= 1;                     % estimate observation model {gamma}

Sim.SuppressGraphics=0;              % suppress graphics output 
Sim.minVar=[0.1 5];                  % min variance limits for kernel dens est

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize barrier and wiener filter parameters
P.rate  = 2;                         % expected spike rate
P.A     = 50;                        % jump size ($\mu$M)
P.tau   = 1.0;                       % calcium decay time constant (sec)
P.lam   = Sim.T/(P.rate*P.A)*Sim.dt; % expected jump size ber time bin
P.sig   = P.A/3;                     % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = -5;                    % linear filter (spike log-rate)
P.C_0       = P.A/2;                 % baseline [Ca++]
P.n         = 1;                     % hill equation exponent
P.k_d       = 100^P.n;               % hill coefficient
P.alpha     = 1;                     % F_max
P.beta      = 0;                     % F_min
P.gamma     = 1e-5;                  % scaled variance
P.zeta      = 4*P.gamma;             % constant variance
P.C_init    = P.C_0;                 % initial [Ca++]
P.tau_c     = P.tau;
P.sigma_c   = P.sig;
P.a         = Sim.dt/P.tau_c;

if Sim.M==1                          % if there are spike history terms
  P.omega = -1;                      % weight
  P.tau_h = 0.02;                    % time constant
  P.sigma_h = 0.01;                  % stan dev of noise
end

P.ignore_lik = 1;                    % ignore josh's-lik-growth cond in EM

P.mismatch_penalty=0.1/Sim.dt;       % mismatch penalty 100ms
P.expected_slack=1.0;                % expected false up to 100%


%% 3) infer nbar and estimate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim.n=double(n>0); Sim.n(Sim.n==0)=NaN;%for plotting purposes
warning off
[I.M I.P]   = GOOPSI_main_v1_0(F,P,Sim);
I.n         = I.M.nbar;
fprintf('\n');
warning on

% evaluate performance over swarm of particles
tgt=cell(1,Sim.N); for i=1:Sim.N tgt{i}=find(n); end;  
est=cell(1,Sim.N); for i=1:Sim.N est{i}=find(I.M.n(i,:)); end
rep=errorsMC(P,est,tgt);             % filter performance report
display(rep)

W=sum(log(max(1e-10,I.M.w)),2);       % show best particle
[a best]=max(W);                


%% 4) prepare summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=[];                                 % summary
c.name  = fname;                      % cell name
c.T     = Sim.T*Sim.dt;               % hold time
c.frame = 1/Sim.dt;                   % frame rate
c.count = sum(n);                     % spikes total
c.rate  = c.count/(Sim.T*Sim.dt);     % spiking rate
c.ntrue = n;                          % true spikes
c.nbar  = I.n;                        % n-bar
c.n     = I.M.n;                      % particles
c.w     = I.M.w;                      % weights
c.P     = I.P;                        % result fit
c.error_dt = rep.position_error * Sim.dt;% spike-time mean error
c.error_fp = rep.false_positives * 100;% percentile false positives
c.error_fn = rep.false_negatives * 100;% percentile false negatives


%% 5) draw sample of spike trains using model with mcmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1=I.P;                               %par structure for SMC
P1.grid=25;                           %size of stochastic grid 
P1.iter=125;                          %number of SMC iterations 

P1.dt=Sim.dt;
P1.sig2=P1.sigma_c^2*P1.dt;
P1.rate=exp(-exp(-P1.k)*P1.dt);

N=size(I.M.C,1);
P1.wf=I.M.w;                          % copy particle swarms     
P1.pf=cell(1,Sim.T);
P1.vf=zeros(2,Sim.T);
for t=1:Sim.T 
  z=[I.M.n(:,t) I.M.C(:,t)]; P1.nf{t}=z;% swarm of particles
  
  v1=0;                               % no variance in n  
  A=repmat(I.M.C(:,t),1,N); A=(A-A'); % [Ca]-[Ca] distances at t
  v2=std(A(:));                       % variance in C as STD of D([Ca]-[Ca])
  P1.vf(:,t)=max(Sim.minVar,[v1 v2])';
end
Y1=cell(1,Sim.T);                     % prep observations
[Y1{:}]=deal(nan);
for t=1:Sim.freq:Sim.T Y1{t}=F(t); end

[Xr Tr]=sampleMCMC(P1,Y1,@spPY_GOOPSI,@spPXX_GOOPSI,@spPG_GOOPSI);%draw sample

K=5;
nf=zeros(floor(P1.iter/K),Sim.T);     % prep bunch of train samples at step K
for k=1:size(nf,1) for t=1:Sim.T nf(k,t)=Tr{K*k}{t}(1); end; end

%evaluate performance over these
m=size(nf,1);
tgt=cell(1,m); for i=1:m tgt{i}=find(n); end;% evaluate filter performance, targets
est=cell(1,m); for i=1:m est{i}=find(nf(i,:)); end% evaluate filter performance, particles
rep=errorsMC(P,est,tgt);             % filter performance report
display(rep)


%% 5) make some diagrams
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rate=exp(P1.k); 
ph=1/P1.gamma;
h2=figure; clf, stem(n,'k'), hold on
stem(0.9*I.M.n(best,:),'c','MarkerSize',8);
stem(0.8*nf(end,:),'m','MarkerSize',8);
plot(I.n,'g','LineWidth',2)
plot(F,'r')
legend('True spikes','Best particle PF','Sample MCMC','nbar PF','Fluorescence','Location','SouthWest')
s=sprintf('Cell R=%gHz ph=%g <dt>=%.3gticks <fp/fn>=%.3g/%.3g',...
  rate,ph,rep.position_error,rep.false_positives,rep.false_negatives);
title(s)

Z=zeros([size(nf) 3]+[6 0 0]);
Z(1,:,1)=n; Z(end,:,1)=n;
Z(2,:,1)=max(0,F)/max(F); Z(end-1,:,1)=max(0,F)/max(F);;
Z(3,:,2)=I.n; Z(end-2,:,2)=I.n;
Z(4:end-3,:,1)=nf;
Z(4:end-3,:,2)=nf;
Z(4:end-3,:,3)=nf;
h3=figure;clf,imagesc(im2uint8(Z)),colormap gray,title(s)
