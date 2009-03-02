% this script fits a cell and prepares parameters for glm-SMC sampler
clear, clc, fprintf('\nfit cell ')
fname=['/Users/joshyv/Research/projects/oopsi/spatial-filter/TS108_6_003.mat']; fprintf('%s\n',fname)
% fname=['../spatial-filter/TS108_6_012.mat']; fprintf('%s\n',fname)
% fname=['../spatial-filter/FilteredData_D080218a24.mat']; fprintf('%s\n',fname)

clear D R
load(fname);

% dt=0.026;
% n=D.n;
% F=0.99*D.F(1:length(n))+0.1;

dt=0.02;                             % time step
F=0.99*D{2}.F+0.01;                  % fluorescence from data
n=zeros(size(F)); n(R.spt)=1;        % spike train from data

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim.T       = length(F);             % # of time steps
Sim.dt      = dt;                    % time step size
Sim.freq    = 1;                     % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;          % # of actual seconds
Sim.T_o     = Sim.T;                 % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;% vector of times
Sim.N       = 100;                   % # of particles
Sim.M       = 0;                     % number of spike history dimensions
Sim.pf      = 1;                     % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                     % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);         % stimulus, for uniform firing

Sim.Mstep   = 1;                     % do M-step
Sim.C_params= 1;                     % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= 1;                     % whether to estimate rate governing parameters {b,k}
Sim.h_params= 1;                     % whether to estimate spike history parameters {h}
Sim.F_params= 1;                     % whether to estimate observation parameters {alpha,beta}
Sim.G_params= 1;                     % whether to estimate observation parameter {gamma}
Sim.MaxIter = 25;                    % # of EM iterartions

Sim.minV=[0 5];                      % min variances in SMC sampler 

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize barrier and wiener filter parameters
P.rate  = 5;                         % expected spike rate
P.A     = 50;                        % jump size ($\mu$M)
P.tau   = 0.5;                       % calcium decay time constant (sec)
P.lam   = Sim.T/(P.rate*P.A)*Sim.dt; % expected jump size ber time bin
P.sig   = P.A/2;                     % standard deviation of noise (\mu M)

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

P.ignorelik = 0;                     % ignore lik-growth cond in EM
P.verbouse  = 1;                     % whether to plot some things

%% 3) infer nbar and estimate parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim.n=double(n>0); Sim.n(Sim.n==0)=NaN;     %for plotting purposes
warning off
[I.M I.P]   = GOOPSI_main_v1_0(F,P,Sim);
I.n         = I.M.nbar;
fprintf('\n');
warning on

tgt=cell(1,100); for i=1:100 tgt{i}=n; end;         % evaluate filter performance, targets
est=cell(1,100); for i=1:100 est{i}=I.M.n(i,:); end % evaluate filter performance, particles
P.mismatch_penalty=0.1/Sim.dt;                      % mismatch penalty 100ms
P.expected_slack=1.0;                               % expected false up to 100%
rep=errorsSMC(P,est,tgt);                           % filter performance report

if(P.verbouse)
  W=sum(log(I.M.w),2); [a b]=max(W);                % show particle with highest weight
  figure,stem(tgt{b}),hold on, stem(est{b}*0.8,'r'), plot(F,'g')
  title(sprintf('Cell %s',fname));
end


%% 4) prepare summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=[];                                               % summary
c.name  = fname;                                    % cell name
c.T     = Sim.T*Sim.dt;                             % hold time
c.frame = 1/Sim.dt;                                 % frame rate
c.count = sum(n);                                   % spikes total
c.rate  = c.count/(Sim.T*Sim.dt);                   % spiking rate
c.nbar  = I.n;                                      % n-bar
c.n     = I.M.n;                                    % particles
c.w     = I.M.w;                                    % weights
c.P     = I.P;                                      % result fit
c.error_dt = rep.time_inaccuracy * Sim.dt;          % spike-time mean error
c.error_fp = rep.false_positives * 100;             % percentile false positives
c.error_fn = rep.false_negatives * 100;             % percentile false negatives
c.error_tgt= tgt{1};                                % match target
c.error_est= est;                                   % match sources

P1=I.P;                                             % par structure for SMC
P1.grid=25;                                         %size of stochastic grid 
P1.iter=250;                                        %number of SMC iterations 
P1.dt=Sim.dt;
P1.sig2=P1.sigma_c^2*P1.dt;
P1.rate=exp(-exp(-P1.k)*P1.dt);
P1.wf=I.M.w;                                        % copy particle swarms     
N=size(I.M.C,1);
P1.pf=cell(1,Sim.T);
P1.vf=zeros(2,Sim.T);
for t=1:Sim.T 
  z=[I.M.n(:,t) I.M.C(:,t)]; P1.nf{t}=z;            % swarm of particles
  
  v1=0;                                             % no variance in n  
  A=repmat(I.M.C(:,t),1,N); A=(A-A');               % [Ca]-[Ca] distances
  v2=std(A(:));                                     % variance in C
  P1.vf(:,t)=max(Sim.minV,[v1 v2])';
end
Y1=cell(1,Sim.T);                                    % convert observations
for t=1:Sim.T Y1{t}=F(t); end



