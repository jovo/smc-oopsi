function [M_best E_best] = GOOPSI_main_v1_0(F,P,Sim)
% this function runs the SMC-EM on a fluorescence time-series, and outputs the inferred
% distributions and parameter estimates
%
% Inputs
% F:    fluorescence time series
% P:    structure of initial parameter estimates
% Sim:  structure of stuff necessary to run smc-em code
%
% Outputs
% M_best:   structure containing mean, variance, and percentiles of inferred distributions
% E_best:   structure containing the final parameter estimates

i           = 0;            % iteration number of EM
k           = 0;            % best iteration so far
P.lik       = -inf;         % we are trying to maximize the likelihood here
maxlik      = P.lik;        % max lik achieved so far
F           = max(F,eps);   % in case there are any zeros in the F time series
Sim.conv    = false;        % EM has NOT yet converged.
Nparticles  = Sim.N;        % store initial Sim parameters

if isfield(Sim,'FastInit')  % if Sim.FastInit=1, then we initialize the parameters using the fast-oopsi code
    if Sim.FastInit==0; FastInit=0; else FastInit=1; end
else
    FastInit=1;
end

if Sim.Mstep==1 && (~isfield(Sim,'SuppressGraphics') || Sim.SuppressGraphics == 0)
    figure(1), clf, nrows=4;
end % if estimating parameters, plot stuff for each iteration

cnt=0;
while Sim.conv==false;
    % some nomenclature to make code easier to read/write
    % these abbrev's are used in forward_step and backward_step
    i           = i+1;                         % index for the iteration of EM
    P.a         = Sim.dt/P.tau_c;
    P.sig2_c    = P.sigma_c^2*Sim.dt;
    P.kx        = P.k'*Sim.x;
    if Sim.M==1
        P.sig2_h    = P.sigma_h.^2*Sim.dt;
        P.g         = 1-Sim.dt/P.tau_h;
    end

    %% forward step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nT = %g steps',Sim.T)
    if isfield(Sim,'TrueSpk') || (FastInit==1 && i==1)                  % force true spikes hack
        fprintf('\nusing provided spike train, skipping forward step\n')
        if isfield(Sim,'TrueSpk')
            S.n     = Sim.TrueSpk;
            S.C     = filter(1,[1 P.a-1],S.n) + P.a*P.C_0;
        else
           [S.n P2] = FOOPSI_v3_05_01(F',P,Sim);
           S.n      = S.n'/max(S.n);
           S.C      = filter(1,[1 -P2.gam],S.n);               % calcium concentration
           S.n(S.n>.2)=1;
           S.n(S.n<=.2)=0;
        end
        if Sim.M>0
            for m=1:Sim.M
                S.h(1,:,m) = filter(1,[1 P.g(m)-1],S.n);
            end
        end
        Sim.N   = 1;
        S.w_f   = 1+0*S.n;
        S.w_b   = S.w_f;
        S.p     = S.w_f;
        M.n_sampl=S.n;
    else
        fprintf('\nforward step:        ')
        Sim.N = Nparticles;
        S = GOOPSI_forward_v1_0(Sim,F,P);
    end;


    %% backward step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nbackward step:       ')
    Z.oney  = ones(Sim.N,1);                    % initialize stuff for speed
    Z.zeroy = zeros(Sim.N);
    Z.C0    = S.C(:,Sim.T);
    Z.C0mat = Z.C0(:,Z.oney)';

    if Sim.C_params==false                          % if not maximizing the calcium parameters, then the backward step is simple
        for t=Sim.T-Sim.freq-1:-1:Sim.freq+1        % actually recurse backwards for each time step
            Z = GOOPSI_backward_v1_0(Sim,S,P,Z,t);
            S.w_b(:,t-1) = Z.w_b;                   % update forward-backward weights
        end
    else                                            % if maximizing calcium parameters,
        % need to compute some sufficient statistics
        M.Q = zeros(3);                             % the quadratic term for the calcium par
        M.L = zeros(3,1);                           % the linear term for the calcium par
        M.J = 0;                                    % remaining terms for calcium par
        M.K = 0;
        for t=Sim.T-Sim.freq-1:-1:Sim.freq+1
            if isfield(Sim,'TrueSpk') || (FastInit==1 && i==1)                  % force true spikes hack
                Z.C0    = S.C(t-1);
                Z.C0mat = Z.C0;
                Z.C1    = S.C(t);
                Z.C1mat = Z.C1;
                Z.PHH   = 1;
                Z.w_b   = 1;
                Z.n1 = S.n(t);
            else
                Z = GOOPSI_backward_v1_0(Sim,S,P,Z,t);
            end
            S.w_b(:,t-1) = Z.w_b;

            % below is code to quickly get sufficient statistics
            C0dt    = Z.C0*Sim.dt;
            bmat    = Z.C1mat-Z.C0mat';
            bPHH    = Z.PHH.*bmat;

            M.Q(1,1)= M.Q(1,1) + sum(Z.PHH*(C0dt.^2));  % Q-term in QP
            M.Q(1,2)= M.Q(1,2) - Z.n1'*Z.PHH*C0dt;
            M.Q(1,3)= M.Q(1,3) + sum(sum(-Z.PHH.*Z.C0mat'*Sim.dt^2));
            M.Q(2,2)= M.Q(2,2) + sum(Z.PHH'*(Z.n1.^2));
            M.Q(2,3)= M.Q(2,3) + sum(sum(Z.PHH(:).*repmat(Z.n1,Sim.N,1))*Sim.dt);
            M.Q(3,3)= M.Q(3,3) + sum(Z.PHH(:))*Sim.dt^2;

            M.L(1)  = M.L(1) + sum(bPHH*C0dt);          % L-term in QP
            M.L(2)  = M.L(2) - sum(bPHH'*Z.n1);
            M.L(3)  = M.L(3) - Sim.dt*sum(bPHH(:));

            M.J     = M.J + sum(Z.PHH(:));              % J-term in QP /sum J^(i,j)_{t,t-1}/

            M.K     = M.K + sum(Z.PHH(:).*bmat(:).^2);  % K-term in QP /sum J^(i,j)_{t,t-1} (d^(i,j)_t)^2/
        end
        M.Q(2,1) = M.Q(1,2);                          % symmetrize Q
        M.Q(3,1) = M.Q(1,3);
        M.Q(3,2) = M.Q(2,3);
    end

    fprintf('\n')

    % copy particle swarm for later
    M.w=S.w_b;
    M.n=S.n;
    if(isfield(S,'h')) M.h=S.h; end
    M.C=S.C;

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

    %% M step

    if Sim.Mstep
        Eold = P;                           % store most recent parameter structure
        P    = GOOPSI_Mstep_v1_0(Sim,S,M,P,F);% update parameters
        fprintf('\n\nIteration #%g, lik=%g, dlik=%g\n',i,P.lik,P.lik-Eold.lik)

        % keep record of best stuff, or if told to ignore lik
        if((isfield(P,'ignorelik') && P.ignorelik==1) || P.lik>= maxlik)
            E_best  = P;                    % update best parameters
            M_best  = M;                    % update best moments
            maxlik  = P.lik;                % update best likelihood
            k       = i;                    % save iteration number of best one
            if(~isfield(Sim,'SuppressGraphics') || ~Sim.SuppressGraphics)
                subplot(nrows,1,4), cla,hold on,% plot spike train estimate
                if isfield(Sim,'n'), stem(Sim.n,'Marker','.',...
                        'MarkerSize',20,'LineWidth',2,'Color',[.75 .75 .75]); end
                
                M.nbar = sum(S.w_b.*S.n,1);
                nvar = sum((repmat(M.nbar,Sim.N,1)-S.n).^2)/Sim.N;
                BarVar=M.nbar+nvar; BarVar(BarVar>1)=1;
                stem(BarVar,'Marker','none','LineWidth',2,'Color',[.8 .8 0]);
                stem(M.nbar,'Marker','none','LineWidth',2,'Color',[0 .5 0])
                axis([0 Sim.T 0 1]),
            end
        end

        % when estimating calcium parameters, display param estimates and lik
        if Sim.C_params==1
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
        if(~isfield(Sim,'SuppressGraphics') || ~Sim.SuppressGraphics)
            if Sim.n_params == true
                fprintf('\nk      = %.2f',P.k)
            end
            subplot(nrows,1,1), hold on, plot(i,P.lik,'o'), axis('tight')
            subplot(nrows,1,2), plot(F,'k'), hold on,
            Cbar = sum(S.w_b.*S.C,1);
            plot(P.alpha*Hill_v1(P,Cbar)+P.beta,'b'), hold off, axis('tight')
            subplot(nrows,1,3), cla, hold on,   % plot spike train estimate
            axis([0 Sim.T 0 1]),
            if isfield(Sim,'n'),
                stem(Sim.n,'Marker','.','MarkerSize',20,'LineWidth',2,...
                    'Color',[.75 .75 .75],'MarkerFaceColor','k','MarkerEdgeColor','k');
                axis('tight'),
            end
            Cvar = sum((repmat(Cbar,Sim.N,1)-S.C).^2)/Sim.N;
            BarVar=M.nbar+nvar; BarVar(BarVar>1)=1;
            stem(BarVar,'Marker','none','LineWidth',2,'Color',[.8 .8 0]);
            stem(M.nbar,'Marker','none','LineWidth',2,'Color',[0 .5 0])

            drawnow
        end

        if i>=Sim.MaxIter
            Sim.conv=true;
        end

    else
        M_best  = M;                     % required for output of function
        E_best  = P;
        Sim.conv= true;
    end

    E_best  = P;                     % update best parameters
    M_best  = M;                     % update best moments

end
fprintf('\n')
