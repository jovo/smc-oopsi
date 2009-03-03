% this script does makes a figure showing how the smc-em approach
% can estimate the linear filter better than other approaches
%
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, fprintf('\nParameter Estimate Simulation Fig\n')

%% 1) set simulation metadata
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Sim.dt      = 1/100;                                % time step size
Sim.freq    = 1;                                    % # of time steps between observations
Sim.N       = 100;                                  % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions

Sim.Mstep   = true;                                 % do M-step
Sim.C_params= true;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params= false;                                % whether to estimate rate governing parameters {b,k}
Sim.h_params= false;                                % whether to estimate spike history parameters {h}
Sim.F_params= false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 49;                                    % max # of EM iterartions
Sim.Plot    = true;                                 % display results for each iteration

%% 2) initialize parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize barrier and wiener filter parameters
% spt     = [108 130 132 136 139 153 171 172 182 194 195 200 222 239 246 252 258 265 273];
P.rate  = 5; %length(spt);                              % expected spike rate
% P.A     = 1;                                        % jump size ($\mu$M)
% P.tau   = 2;                                        % calcium decay time constant (sec)
% P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-P.rate*Sim.dt)/Sim.dt);    % linear filter
% P.k         = [1 -P.k P.k -2*P.k 2*P.k]';
P.tau_c     = 0.5;
P.A         = 5;
P.C_0       = 5;                                   % baseline [Ca++]
P.C_init    = P.C_0;                                 % initial [Ca++]
P.sigma_c   = 1;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.gamma     = .3e-5;                                 % scaled variance
P.zeta      = 4*P.gamma;                            % constant variance
P.a         = Sim.dt/P.tau_c;

if Sim.M==1                                         % if there are spike history terms
    P.omega = -1;                                   % weight
    P.tau_h = 0.5;                                  % time constant
    P.sigma_h = 0.01;                               % stan dev of noise
end

%% 3) simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Algs    = 9;%[4 9];                                    % which algorithms within DataComp to use
Sim.T   = 5000;%[100 500 3000 6000];
Nspikes = [5 10 20 40];
fig     = figure(100); clf,
wh      = [7 5]; set(fig,'PaperPosition',[0 11-wh(2) wh]);
nrows   = length(Nspikes);
gray    = [.75 .75 .75];                            % define gray color
col     = [1 0 0; 0 .5 0; 0 0 1; 1 0 1];            % define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;                   % define colors for std
inter   = 'tex';                                    % interpreter for axis labels
fs      = 12;                                       % font size
ms      = 20;                                       % marker size for real spike
sw      = 3.5;                                      % spike width
lw      = 2;                                        % line width
j       = 100;                                      % ln(j) is the # of sig digs to display in legend

for T=1:length(Nspikes)
    for tr=1:10                                     % do a number of trials
        fprintf('\nTrial %.0f\n',tr)
        Sim.T       = 5000;
        Sim.T_o     = Sim.T;                        % # of observations
        tvec        = Sim.freq:Sim.freq:Sim.T;      % vector of observation times
        Sim.x       = 1*randn(Sim.StimDim,Sim.T);   % stimulus
        Sim.x(1,:)  = 1;                            % for baseline firing rate
        P.lam       = Sim.T/(P.rate*P.A)*Sim.dt;    % expected jump size ber time bin
        kx          = P.k'*Sim.x;                   % external input to neuron
        U_sampl     = rand(1,Sim.T);                % generate random number to use for sampling
        if Sim.M==1                                 % if spike history terms, recursively
            p         = zeros(1,Sim.T);             % extize p_t because it must be updated iteratively
            n         = zeros(1,Sim.T);             % extize n_t because it must be updated iteratively
            h         = zeros(1,Sim.T);             % extize spike history because it must be updated iteratively
            epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(1,Sim.T); % generate noise on spike history
            for t=2:Sim.T                           % update states
                h(:,t)  = (1-Sim.dt./P.tau_h).*h(:,t-1)+n(t-1) + epsilon_h(:,t);% update h terms
                y_t     = kx(t)+P.omega'*h(:,t);    % generate operand for rate function
                p(t)    = 1-exp(-exp(y_t)*Sim.dt);  % generate rate
                n(t)    = U_sampl(t)<p(t);          % sample from bernoulli with prob p_t
            end %time loop
        else
            p=1-exp(-exp(P.k'*Sim.x)*Sim.dt);
            n=(rand(Sim.T,1)<p')';
        end

        spt = find(n,Nspikes(T)+1);         % store spike train
        Sim.T = max(spt)-1;                 % cut time down to include only a certain number of spikes
        Sim.x = 1*randn(Sim.StimDim,Sim.T);         % stimulus
        Sim.x(1,:)  = 1;                            % for baseline firing rate

        C=zeros(Sim.T,1);
        epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
        for t=2:Sim.T                                           %recursively update calcium
            C(t)  = (1-Sim.dt/P.tau_c)*C(t-1) + P.a*P.C_0+ P.A*n(t) + epsilon_c(t);
        end
        % C=filter(P.A,[1 -(1-Sim.dt/P.tau)],n)+P.C_0;
        S=Hill_v1(P,C)';
        F=(P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(1,Sim.T))';
        F(F<=0)=eps;
        figure(100*T+tr), clf, plot(z1(F)+1), hold on, stem(n(1:Sim.T));
        Sim.n = double(n(1:Sim.T)); Sim.n(Sim.n==0)=NaN;     % for plotting purposes in ParticleFiltD

        D{tr,T}.F   = F;                            % store fluorescence
        D{tr,T}.spt = spt(1:end-1);
        %         D{tr,T}.k   = (Sim.x*n')/sum(n);            % k estimated using actual spike trains
        %         D{tr,T}.Fk  = (Sim.x(:,tvec)*F(tvec))/sum(F(tvec));% k estimated using raw fluorescence signal

        %% 4) infer spikes and estimate parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for m=Algs
            fprintf('\nAlg %X\n',m)
            E = P; if tr<=5, aa=.5; else aa=2; end
            if Sim.Mstep == true;
                E.sigma_c=10*E.sigma_c;
                if Sim.C_params == true, E.A = E.A*aa; E.tau_c = E.tau_c*aa; E.C_0 = E.C_0*aa; E.sig = 20; end
                if Sim.n_params == true, E.k = E.k./(aa*ones(size(E.k))); end
                if Sim.h_params == true, E.omega = E.omega*aa; end
                if Sim.F_params == true, E.alpha = E.alpha*aa; E.beta = E.beta*aa; end
            end
            Sim.Alg = m;
            if m==4
                G=F(tvec);
                Tim=Sim; Tim.dt=Tim.dt*Tim.freq;
                I{m,tr,T}  = DataComp13(G,E,Tim);
                I{m,tr,T}.P.k=(Tim.x(:,tvec)*I{m,tr,T}.n)/sum(I{m,tr,T}.n);
            else
                I{m,tr,T}  = DataComp13(F,E,Sim);
            end
        end
        save
    end % trials


    %% 5) plot results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    figure(T), clf
    for tr=1:10
        %     for T=1:3
        A(tr,T)=I{m,tr,T}.P.A;
        tau(tr,T)=I{m,tr,T}.P.tau_c;
        C_0(tr,T)=I{m,tr,T}.P.C_0;
        D{tr,T}.n=zeros(size(D{tr,T}.F));
        D{tr,T}.n(D{tr,T}.spt)=1;
        subplot(10,1,tr),
        plot(z1(D{tr,T}.F),'k'), hold on, stem(D{tr,T}.n), bar(I{m,tr,T}.M.nbar), axis('tight')
        %     end
    end
    drawnow
    
end % Time steps
set(gca,'XTick',0:6,'XTickLabel',0:6,'YTick',[-1 0 1])