% this script fits a cell and prepares parameters for glm-SMC sampler
clear, close all, clc, fprintf('\nfit cell ')

% dt=0.026;

T=10;
dt=0.02;                             % time step

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sim.T       = round(T/dt);           % # of time steps
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
Sim.MaxIter = 1;                    % # of EM iterartions

Sim.minV=[0.1 2.5];                  % min variances in SMC sampler 

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SNR=2; SNR_BG=5; ph=1e5;
% initialize barrier and wiener filter parameters
P.rate  = 2.5;                       % expected spike rate
P.A     = 25;                        % jump size ($\mu$M)
P.tau   = 1.0;                       % calcium decay time constant (sec)
P.lam   = Sim.T/(P.rate*P.A)*Sim.dt; % expected jump size ber time bin
P.sig   = P.A/SNR;                   % standard deviation of noise (\mu M)

% initialize particle filter parameters
rf=log(-log(1-P.rate/Sim.T)/Sim.dt);
P.k         = rf;                    % linear filter (spike log-rate)
P.C_0       = P.A/SNR_BG;            % baseline [Ca++]
P.n         = 1;                     % hill equation exponent
P.k_d       = 100^P.n;               % hill coefficient
P.alpha     = 1;                     % F_max
P.beta      = 0;                     % F_min
P.gamma     = 1/ph;                  % scaled variance
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

%% 3) simulate cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n=rand(1,Sim.T)>exp(-P.rate*Sim.dt);
C=zeros(1,Sim.T); C(1)=P.C_0;
for t=2:Sim.T C(t)=max(0,C(t-1)+Sim.dt/P.tau*(P.C_0-C(t-1))+P.A*n(t)+P.sig*sqrt(Sim.dt)*randn); end
S=Hill_v1(P,C);
F=P.beta+P.alpha*S+sqrt(P.zeta+P.gamma*S).*randn(size(S));

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
display(rep);

if(P.verbouse)
  W=sum(log(I.M.w),2); [a best]=max(W);                % show particle with highest weight
  figure,stem(tgt{best}),hold on, stem(est{best}*0.8,'r'), plot(F,'g')
end


%% 4) prepare summary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P1=I.P;                                             % par structure for SMC
P1.grid=25;                                         % size of stochastic grid 
P1.iter=10;                                         % number of SMC iterations 
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

[Xr Tr]=sampleSMC(P1,Y1,@spPY_GOOPSI,@spPXX_GOOPSI,@spPG_GOOPSI);%draw sample

nf=zeros(floor(P1.iter/5),Sim.T);
for k=1:size(nf,1) for t=1:Sim.T nf(k,t)=Tr{5*k}{t}(1); end; end

%% evaluate errors in this
nsamples=2;
tgt=cell(1,nsamples); for i=1:nsamples tgt{i}=n; end;         % evaluate filter performance, targets
est=cell(1,nsamples); for i=1:nsamples est{i}=nf(i,:); end    % evaluate filter performance, particles
P.mismatch_penalty=0.1/Sim.dt;                      % mismatch penalty 100ms
P.expected_slack=1.0;                               % expected false up to 100%
rep=errorsSMC(P,est,tgt);                           % filter performance report
display(rep)

%make some graphs
figure,stem(n,'k'),hold on
plot(F,'r')
plot(I.n)
stem(0.9*I.M.n(best,:));
stem(0.8*nf(end,:));
legend('True spikes','F','nbar','Best PF train','SMC train')

Z=zeros([size(nf) 3]);
Z(:,:,1)=repmat(n,size(nf,1),1);
Z(:,:,2)=nf;
figure,imagesc(nf),title('SMC trains raster')



