function [S M] = smc_em_bern_FoBaMo_v5(Sim,R,P)
% this is really the main best of the code --- it does the E step for our
% SMC-EM algorithm.  

%% makes code prettier: this just makes the code easier to read and write

B.sig2_c    = P.sigma_c^2*Sim.dt;
B.a         = Sim.dt/P.tau_c;
B.A         = P.A;
B.kx        = P.k'*Sim.x;
if Sim.M>0
    B.sig2_h    = P.sigma_h.^2*Sim.dt;
    B.g         = 1-Sim.dt/P.tau_h;
    B.omega     = P.omega;
end
B.C_0       = P.C_0;
B.C_init    = P.C_init;
B.n         = P.n;
B.k_d       = P.k_d;
B.alpha     = P.alpha;
B.beta      = P.beta;
B.gamma     = P.gamma;
B.zeta      = P.zeta;

%% forward step
fprintf('\nTotal Number of Steps=%d\n',Sim.T)   %print to screen the total number of time steps in the simulation
if Sim.pf==0                                    %if Sim.pf==0, then use the PRIOR sampler
    fprintf('\nprior forward step:          ')
    S = smc_em_bern_PriorSampl_v4(Sim,R,B);
elseif Sim.pf==1                                %else, use the conditional sampler
    fprintf('\nk-mixture forward step:      ')
    S = smc_em_bern_k_comp_BackSampl_v3(Sim,R,B);
end

%% backward step
fprintf('\nbackward step:               ')
Z.oney  = ones(Sim.N,1);                        %initialize stuff for speed
Z.zeroy = zeros(Sim.N);
Z.C0    = S.C(:,Sim.T);
Z.C0mat = Z.C0(:,Z.oney)';

if Sim.Mstep==false                             %if not maximizing the parameters, then the backward step is simple
    for t=Sim.T-Sim.freq-1:-1:Sim.freq+1        %actually recurse backwards for each time step
        Z = smc_em_bern_backwardPF_v10(Sim,S,B,Z,t);
        S.w_b(:,t-1) = Z.w_b;                   %update forward-backward weights
    end
else                                            %if maximizing the parameters, we need to also compute some other sufficient statistics
    M.Q = zeros(3);                             %the quadratic term for the calcium parameters
    M.L = zeros(3,1);                           %the linear term for the calcium parameters
    M.u = 0;                                    %not used
    for t=Sim.T-Sim.freq-1:-1:Sim.freq+1
        Z = smc_em_bern_backwardPF_v10(Sim,S,B,Z,t);
        S.w_b(:,t-1) = Z.w_b;

        % below is code to quickly and recursively get the suff. stats.
        % necessary to estimate the calcium parameters
        C0dt    = Z.C0*Sim.dt;
        bmat    = Z.C1mat-Z.C0mat';
        bPHH    = Z.PHH.*bmat;
        syms A a C_0

        M.Q(1,1)= M.Q(1,1) + sum(Z.PHH*(C0dt.^2));
        M.Q(1,2)= M.Q(1,2) - Z.n1'*Z.PHH*C0dt;
        M.Q(1,3)= M.Q(1,3) + sum(sum(-Z.PHH.*Z.C0mat'*Sim.dt^2));
        M.Q(2,2)= M.Q(2,2) + sum(Z.PHH'*(Z.n1.^2));
        M.Q(2,3)= M.Q(2,3) + sum(sum(Z.PHH(:).*repmat(Z.n1,Sim.N,1))*Sim.dt);
        M.Q(3,3)= M.Q(3,3) + sum(Z.PHH(:))*Sim.dt^2;

        M.L(1)  = M.L(1) + sum(bPHH*C0dt);
        M.L(2)  = M.L(2) - sum(bPHH'*Z.n1);
        M.L(3)  = M.L(3) - Sim.dt*sum(bPHH(:));
    end
    M.Q(2,1) = M.Q(1,2);
    M.Q(3,1) = M.Q(1,3);
    M.Q(3,2) = M.Q(2,3);
end

fprintf('\n')

%%   compute moments
M.nbar = sum(S.w_b.*S.n,1);
M.nvar = sum((repmat(M.nbar,Sim.N,1)-S.n).^2)/Sim.N;

M.Cbar = sum(S.w_b.*S.C,1);
M.Cvar = sum((repmat(M.Cbar,Sim.N,1)-S.C).^2)/Sim.N;

%%   compute mse
% M.mse.n=sum(S.w_b.*(S.n-repmat(R.n,Sim.N,1)).^2);
% M.mse.C=sum(S.w_b.*(S.C-repmat(R.C,Sim.N,1)).^2);
