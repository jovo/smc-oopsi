%this script runs simulations and prepares some plots
clear, close all, clc, fprintf('\nfit cell ')

gallery=[];
rate=10;                    % assumed firing rate
SNR=3.5;                    % assumed A_[Ca]/Noise[Ca]
BGR=4;                      % assumed A_[Ca]/BGR[Ca]

for frame_rate=[100 60 30 10]% do these frame rates
  for ph=[1e5 1e4 1e3]       % do these photo counts

    if(rate  > frame_rate) continue; end%don't do this

    T=max(5,25/rate);       %recording duration 5s but at least 25 spikes
    dt=1/frame_rate;        % frame time step

    %% 1) set simulation metadata
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Sim.T       = round(T/dt);           % # of time steps
    Sim.dt      = dt;                    % time step size
    Sim.freq    = 1;                     % # of time steps between observations
    Sim.Nsec    = Sim.T*Sim.dt;          % # of actual seconds
    Sim.T_o     = Sim.T;                 % # of observations
    Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;% vector of times
    Sim.N       = 50;                    % # of particles
    Sim.M       = 0;                     % number of spike history dimensions
    Sim.pf      = 1;                     % use conditional sampler (not prior) when possible
    Sim.StimDim = 1;                     % # of stimulus dimensions
    Sim.x       = ones(1,Sim.T);         % stimulus, for uniform firing

    Sim.Mstep   = 1;                     % do M-step
    Sim.MaxIter = 15;                    % # of EM iterartions
    Sim.C_params= 1;                     % estimate calcium model {tau,A,C_0,sig}
    Sim.n_params= 1;                     % estimate rate model {b,k}
    Sim.h_params= 1;                     % estimate spike history model {h}
    Sim.F_params= 1;                     % estimate observation model {alpha,beta}
    Sim.G_params= 1;                     % estimate observation model {gamma}
    
    Sim.Scan    = 0;                     % scans or epi data
    Sim.SuppressGraphics=0;              % suppress graphics output
    Sim.minVar=[0.1 5];                  % min variance limits for kernel dens est

    %% 2) initialize parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('Cell [%g %g %g %g]\n',SNR,frame_rate,ph,rate);

    % initialize barrier and wiener filter parameters
    P.rate  = rate;                       % expected spike rate
    P.A     = 35;                         % jump size ($\mu$M)
    P.tau   = 0.25;                       % calcium decay time constant (sec)
    P.lam   = Sim.T/(P.rate*P.A)*Sim.dt;  % expected jump size ber time bin
    P.sig   = P.A/SNR;                    % standard deviation of noise (\mu M)

    % initialize particle filter parameters
    rf          = log(-log(1-P.rate/Sim.T)/Sim.dt);
    P.k         = rf;                    % linear filter (spike log-rate)
    P.C_0       = P.A/BGR;               % baseline [Ca++]
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
    P.mismatch_penalty=0.1/Sim.dt;       % mismatch penalty 100ms
    P.expected_slack=1.0;                % expected false up to 100%

    %% 3) simulate cell
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ntgt=[];    
    n=poissrnd(P.rate*Sim.dt,1,Sim.T);   % POISSON GENERATION MODEL
    for t=find(n) ntgt=[ntgt,repmat(t,1,n(t))]; end% clone multiple spikes
    if(length(ntgt)==0) m=m-1; continue; end% try again -- no spikes
    
    C=zeros(1,Sim.T); 
    C(1)=max(0,P.C_0+P.sig*randn);
    for t=2:Sim.T                        % setup due Vogelstein-Paninski
      C(t)=max(0,C(t-1)+Sim.dt/P.tau*(P.C_0-C(t-1))+P.A*n(t)+P.sig*sqrt(Sim.dt)*randn);
    end
    S=Hill_v1(P,C);                      % fluorescence
    F=P.beta+P.alpha*S+sqrt(P.zeta+P.gamma*S).*randn(size(S));

    %% 3) infer nbar and estimate parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initialize barrier and wiener filter parameters - MESS UP PAR AS IF WE DON'T KNOW
    P.rate  = 1.0;                       % expected spike rate
    P.A     = 50;                        % jump size ($\mu$M)
    P.tau   = 0.5;                       % calcium decay time constant (sec)
    P.lam   = Sim.T/(P.rate*P.A)*Sim.dt; % expected jump size ber time bin
    P.sig   = P.A/3;                     % standard deviation of noise (\mu M)
    P.C_0       = P.A/3;                 % baseline [Ca++]
    P.alpha     = 1;                     % F_max
    P.beta      = 0;                     % F_min
    P.gamma     = 1e-5;                  % scaled variance
    P.zeta      = 4*P.gamma;             % constant variance
    P.C_init    = P.C_0;                 % initial [Ca++]
    P.tau_c     = P.tau;
    P.sigma_c   = P.sig;
    P.a         = Sim.dt/P.tau_c;

    nz=zeros(size(F)); nz(ntgt)=1;
    Sim.n=double(nz>0); Sim.n(Sim.n==0)=NaN;%for plotting purposes
    warning off
    [I.M I.P]   = GOOPSI_main_v1_0(F,P,Sim);
    I.n         = I.M.nbar;
    fprintf('\n');
    warning on

    W=sum(log(max(1e-100,I.M.w)),2); 
    [a best]=max(W);   % best particle    
    
    %% 4) prepare summary
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    P1=I.P;                              % par structure for SMC
    P1.grid=25;                          % size of stochastic grid
    P1.iter=125;                         % number of SMC iterations
    
    P1.dt=Sim.dt;
    P1.sig2=P1.sigma_c^2*P1.dt;
    P1.rate=exp(-exp(-P1.k)*P1.dt);
    
    P1.wf=I.M.w;                         % copy particle swarms
    N=size(I.M.C,1);
    P1.pf=cell(1,Sim.T);
    P1.vf=zeros(2,Sim.T);
    for t=1:Sim.T
      z=[I.M.n(:,t) I.M.C(:,t)]; P1.nf{t}=z;% swarm of particles

      v1=0;                                 % no variance in n
      A=repmat(I.M.C(:,t),1,N); A=(A-A');   % [Ca]-[Ca] distances
      v2=max(P1.sigma_c,std(A(:)));         % variance as Var(D[Ca]-[Ca])
      P1.vf(:,t)=max(Sim.minVar,[v1 v2])';
    end
    Y1=cell(1,Sim.T);                       % convert observations
    [Y1{:}]=deal(nan);                      % skipped observations
    for t=1:Sim.freq:Sim.T Y1{t}=F(t); end

    %MCMC
    [Xr Tr]=sampleMCMC(P1,Y1,@spPY_GOOPSI,@spPXX_GOOPSI,@spPG_GOOPSI);

    K=5;
    nf=zeros(floor(P1.iter/K),Sim.T);       % draw bunch of samples every K steps
    for k=1:size(nf,1) for t=1:Sim.T nf(k,t)=Tr{K*k}{t}(1); end; end

    %evaluate errors in this
    mf=size(nf,1);                          % evaluate performance, targets and samples
    tgt=cell(1,mf); for i=1:mf tgt{i}=ntgt(ntgt>5 & ntgt+5<Sim.T); end% avoid border spk
    est=cell(1,mf); for i=1:mf est{i}=find(nf(i,:)); end
    for i=1:mf est{i}=est{i}(est{i}>5 & est{i}+5<Sim.T); end
    rep=errorsMC(P,est,tgt);               % performance report
    display(rep)

    %make some graphs
    if(Sim.SuppressGraphics == 0)
      h2=figure(2);clf,stem(nz,'k'),hold on
      stem(0.9*I.M.n(best,:),'c','MarkerSize',8);
      stem(0.8*nf(end,:),'m','MarkerSize',8);
      plot(I.n,'g','LineWidth',2)
      plot(F,'r')
      legend('True spikes','Best particle PF','Sample MCMC','nbar PF','Fluorescence','Location','SouthWest')
      s=sprintf('Cell R=%gHz FR=%gHz F=%iK <dt>=%.3gticks <fp>/<fn>=%.3g/%.3g',...
        rate,frame_rate,round(ph/1000),rep.position_error*P1.dt,rep.false_positives,rep.false_negatives);
      title(s)

      Z=zeros([size(nf) 3]+[6 0 0]);
      Z(1,:,1)=nz; Z(end,:,1)=nz;
      Z(2,:,1)=max(0,F)/max(F); Z(end-1,:,1)=max(0,F)/max(F);;
      Z(3,:,2)=I.n; Z(end-2,:,2)=I.n;
      Z(4:end-3,:,1)=nf;
      Z(4:end-3,:,2)=nf;
      Z(4:end-3,:,3)=nf;
      h3=figure(3);clf,imagesc(im2uint8(Z)),colormap gray,title(s)

      %SAVE FILE HERE
      fname=sprintf('simCOSHIGH%gHz%gHz%iK',rate,frame_rate,round(ph/1000));
      hgsave(h2,[fname,'-plot.fig']);
      hgsave(h3,[fname,'-raster.fig']);
      pause(1)
    end
  end
end