function [M_best E_best V] = smc_oopsi(F,V,P)
% this function runs the SMC-EM on a fluorescence time-series, and outputs the inferred
% distributions and parameter estimates
%
% Inputs
% F:    fluorescence time series
% V:    structure of stuff necessary to run smc-em code
% P:    structure of initial parameter estimates
%
% Outputs
% M_best:   structure containing mean, variance, and percentiles of inferred distributions
% E_best:   structure containing the final parameter estimates
% V:        structure Variables for algorithm to run

if nargin < 2,          V       = struct;       end
if ~isfield(V,'T'),     V.T     = length(F);    end     % # of observations
if ~isfield(V,'freq'),  V.freq  = 1;            end     % # time steps between observations
if ~isfield(V,'T_o'),   V.T_o   = V.T;          end     % # of observations
if ~isfield(V,'N'),     V.N     = 99;           end     % # particles
if ~isfield(V,'M'),     V.M     = 0;            end     % # of spike history terms
if ~isfield(V,'pf'),    V.pf    = 1;            end     % whether to use conditional sampler
if ~isfield(V,'x'),     V.x     = ones(1,V.T);  end     % stimulus
if ~isfield(V,'Scan'),  V.Scan  = 0;            end     % epi or scan
if ~isfield(V,'name'),  V.name  ='oopsi';       end     % name for output and figure
if ~isfield(V,'ignorelike'),  V.ignorelik  = 1; end     % epi or scan
if ~isfield(V,'true_n'),V.true_n = 0;           end     % whether true spikes are available   
if ~isfield(V,'smc_iter_max'),                          % max # of iterations before convergence
    reply = str2double(input('\nhow many EM iterations would you like to perform \nto estimate parameters (0 means use default parameters): ', 's'));
    V.smc_iter_max = reply;
end
if ~isfield(V,'dt'),
    fr = input('what was the frame rate for this movie (in Hz)? ');
    V.dt = 1/fr;
end

% set which parameters to estimate
if ~isfield(V,'est_c'),     V.est_c   = 1;      end     % calcium params
if ~isfield(V,'est_n'),     V.est_n   = 1;      end     % b,k
if ~isfield(V,'est_h'),     V.est_h   = 0;      end     % w
if ~isfield(V,'est_F'),     V.est_F   = 1;      end     % alpha, beta
if ~isfield(V,'smc_plot'),  V.smc_plot = 1;     end     % plot results with each iteration
if ~isfield(V,'holdTau'),   V.holdTau=0;        end     % whether to hold tau_c fixed while estimating parameters (useful when data is poor)
if isfield(V,'true_n'),     true_n=1;           else true_n=0;  end
if V.smc_iter_max>1 && V.smc_plot == 1                  % if estimating parameters, plot stuff for each iteration
    figNum=10;
    figure(figNum), clf, nrows=4;
end

%% initialize model Parameters

if nargin < 3,          P       = struct;       end
if ~isfield(P,'tau_c'), P.tau_c = 1;            end     % calcium decay time constant (sec)
if ~isfield(P,'A'),     P.A     = 50;           end     % change ins [Ca++] after a spike (\mu M)
if ~isfield(P,'C_0'),   P.C_0   = 0;            end     % baseline [Ca++] (\mu M)
if ~isfield(P,'C_init'),P.C_init= 0;            end     % initial [Ca++] (\mu M)
if ~isfield(P,'sigma_c'),P.sigma_c= 0.1;        end     % standard deviation of noise (\mu M)
if ~isfield(P,'n'),     P.n     = 1;            end     % hill equation exponent
if ~isfield(P,'k_d'),   P.k_d   = 200;          end     % hill coefficient
if ~isfield(P,'k'),                                     % linear filter
    k   = str2double(input('approx. how many spikes underly this trace: ', 's'));
    P.k = log(-log(1-k/V.T)/V.dt); 
end  
if ~isfield(P,'alpha'), P.alpha = mean(F);      end     % scale of F
if ~isfield(P,'beta'),  P.beta  = min(F);       end     % offset of F
if ~isfield(P,'zeta'),  P.zeta  = P.alpha/5;    end     % constant variance
if ~isfield(P,'gamma'), P.gamma = P.zeta/5;     end     % scaled variance
if V.M==1                                               % if there are spike history terms
    if ~isfield(P,'omega'),   P.omega   = -1;   end     % weight
    if ~isfield(P,'tau_h'),   P.tau_h   = 0.02; end     % time constant
    if ~isfield(P,'sigma_h'), P.sigma_h = 0;    end     % stan dev of noise
end

%% initialize other stuff
i           = 0;            % iteration number of EM
i_best      = 0;            % best iteration so far
P.lik       = -inf;         % we are trying to maximize the likelihood here
maxlik      = P.lik;        % max lik achieved so far
F           = max(F,eps);   % in case there are any zeros in the F time series
conv        = false;        % EM has NOT yet converged.
Nparticles  = V.N;          % store initial V parameters
cnt         = 0;
starttime   = cputime;

while conv==false;
    i           = i+1;                                  % index for the iteration of EM
    P.a         = V.dt/P.tau_c;                         % for brevity
    V.smc_iter_tot = i;                                 % total number of iterations completed

    %% forward step
    S   = smc_oopsi_forward(V,F,P);

    %% backward step
    fprintf('\nbackward step:       ')
    Z.oney  = ones(V.N,1);                          % initialize stuff for speed
    Z.zeroy = zeros(V.N);
    Z.C0    = S.C(:,V.T);
    Z.C0mat = Z.C0(:,Z.oney)';

    if V.est_c==false                               % if not maximizing the calcium parameters, then the backward step is simple
        if true_n==1                                % when spike train is provided, backwards is not necessary
            S.w_b=S.w_f;
        else
            for t=V.T-V.freq-1:-1:V.freq+1          % actually recurse backwards for each time step
                Z = smc_oopsi_backward(V,S,P,Z,t);
                S.w_b(:,t-1) = Z.w_b;               % update forward-backward weights
            end
        end
    else                                            % if maximizing calcium parameters,
                                                    % need to compute some sufficient statistics
        M.Q = zeros(3);                             % the quadratic term for the calcium par
        M.L = zeros(3,1);                           % the linear term for the calcium par
        M.J = 0;                                    % remaining terms for calcium par
        M.K = 0;
        for t=V.T-V.freq-1:-1:V.freq+1
            if true_n                               % force true spikes hack
                Z.C0    = S.C(t-1);
                Z.C0mat = Z.C0;
                Z.C1    = S.C(t);
                Z.C1mat = Z.C1;
                Z.PHH   = 1;
                Z.w_b   = 1;
                Z.n1 = S.n(t);
            else
                Z = smc_oopsi_backward(V,S,P,Z,t);
            end
            S.w_b(:,t-1) = Z.w_b;

            % below is code to quickly get sufficient statistics
            C0dt    = Z.C0*V.dt;
            bmat    = Z.C1mat-Z.C0mat';
            bPHH    = Z.PHH.*bmat;

            M.Q(1,1)= M.Q(1,1) + sum(Z.PHH*(C0dt.^2));  % Q-term in QP
            M.Q(1,2)= M.Q(1,2) - Z.n1'*Z.PHH*C0dt;
            M.Q(1,3)= M.Q(1,3) + sum(sum(-Z.PHH.*Z.C0mat'*V.dt^2));
            M.Q(2,2)= M.Q(2,2) + sum(Z.PHH'*(Z.n1.^2));
            M.Q(2,3)= M.Q(2,3) + sum(sum(Z.PHH(:).*repmat(Z.n1,V.N,1))*V.dt);
            M.Q(3,3)= M.Q(3,3) + sum(Z.PHH(:))*V.dt^2;

            M.L(1)  = M.L(1) + sum(bPHH*C0dt);          % L-term in QP
            M.L(2)  = M.L(2) - sum(bPHH'*Z.n1);
            M.L(3)  = M.L(3) - V.dt*sum(bPHH(:));

            M.J     = M.J + sum(Z.PHH(:));              % J-term in QP /sum J^(i,j)_{t,t-1}/

            M.K     = M.K + sum(Z.PHH(:).*bmat(:).^2);  % K-term in QP /sum J^(i,j)_{t,t-1} (d^(i,j)_t)^2/
        end
        M.Q(2,1) = M.Q(1,2);                            % symmetrize Q
        M.Q(3,1) = M.Q(1,3);
        M.Q(3,2) = M.Q(2,3);
    end
    fprintf('\n')

    % copy particle swarm for later
    M.w = S.w_b;
    M.n = S.n;
    M.C = S.C;
    if isfield(S,'h'), M.h=S.h; end

    % check failure mode caused by too high P.A (low P.sigma_c)
    fact=1.55;
    if(sum(S.n(:))==0 && cnt<10)                % means no spikes anywhere
        fprintf(['Failed to find any spikes, likely too high a P.A.\n',...
            'Attempting to lower by factor %g...\n'],fact);
        P.A=P.A/fact;
        P.C_0=P.C_0/fact;
        P.sigma_c=P.sigma_c/fact;
        cnt=cnt+1;
        continue;
    elseif(cnt>=10)
        M_best=M;
        E_best=P;
        fprintf('Warning: there are no spikes in the data. Wrong initialization?');
        return;
    end
    M.nbar = sum(S.w_b.*S.n,1);

    %% M step
    if V.smc_iter_max>1
        Eold = P;                           % store most recent parameter structure
        P    = smc_oopsi_m_step(V,S,M,P,F); % update parameters
        fprintf('\n\nIteration #%g, lik=%g, dlik=%g\n',i,P.lik,P.lik-Eold.lik)

        % keep record of best stuff, or if told to ignore lik
        if V.ignorelik==1 || P.lik>= maxlik
            E_best  = P;                    % update best parameters
            M_best  = M;                    % update best moments
            maxlik  = P.lik;                % update best likelihood
            i_best  = i;                    % save iteration number of best one
            if V.smc_plot
                figure(figNum)
                subplot(nrows,1,4), cla,hold on,% plot spike train estimate
                if isfield(V,'n'), stem(V.n,'Marker','.',...
                        'MarkerSize',20,'LineWidth',2,'Color',[.75 .75 .75]); end
                nvar = sum((repmat(M.nbar,V.N,1)-S.n).^2)/V.N;
                BarVar=M.nbar+nvar; BarVar(BarVar>1)=1;
                stem(BarVar,'Marker','none','LineWidth',2,'Color',[.8 .8 0]);
                stem(M.nbar,'Marker','none','LineWidth',2,'Color',[0 .5 0])
                axis([0 V.T 0 1]),
            end
        end

        % when estimating calcium parameters, display param estimates and lik
        if V.est_c==1
            dtheta  = norm([P.tau_c; P.A; P.C_0]-...
                [Eold.tau_c; Eold.A; Eold.C_0])/norm([Eold.tau_c; Eold.A; Eold.C_0; P.sigma_c]);
            fprintf('\ndtheta = %.2f',dtheta);
            fprintf('\ntau    = %.2f',P.tau_c)
            fprintf('\nA      = %.2f',P.A)
            fprintf('\nC_0    = %.2f',P.C_0)
            fprintf('\nsig    = %.2f',P.sigma_c)
            fprintf('\nalpha  = %.2f',P.alpha)
            fprintf('\nbeta   = %.2f',P.beta)
            fprintf('\ngamma  = %.2g',P.gamma)
        end

        % plot lik and inferrence
        if V.smc_plot
            if V.est_n == true
                fprintf('\nk      = %.2f',P.k)
            end
            subplot(nrows,1,1), hold on, plot(i,P.lik,'o'), axis('tight')
            subplot(nrows,1,2), plot(F,'k'), hold on,
            Cbar = sum(S.w_b.*S.C,1);
            plot(P.alpha*Hill_v1(P,Cbar)+P.beta,'b'), hold off, axis('tight')
            subplot(nrows,1,3), cla, hold on,   % plot spike train estimate
            axis([0 V.T 0 1]),
            if isfield(V,'n'),
                stem(V.n,'Marker','.','MarkerSize',20,'LineWidth',2,...
                    'Color',[.75 .75 .75],'MarkerFaceColor','k','MarkerEdgeColor','k');
                axis('tight'),
            end
            Cvar = sum((repmat(Cbar,V.N,1)-S.C).^2)/V.N;
            BarVar=M.nbar+nvar; BarVar(BarVar>1)=1;
            stem(BarVar,'Marker','none','LineWidth',2,'Color',[.8 .8 0]);
            stem(M.nbar,'Marker','none','LineWidth',2,'Color',[0 .5 0])

            drawnow
        end

        if i>=V.smc_iter_max
            conv=true;
        end

    else
        M_best  = M;                     % required for output of function
        E_best  = P;
        conv    = true;
    end

    %     M_best  = M;                     % update best moments
    %     E_best  = P;                     % update best parameters

end
fprintf('\n')
V.smc_time=cputime-starttime;
V=orderfields(V);