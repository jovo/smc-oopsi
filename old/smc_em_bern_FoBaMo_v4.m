function [S M] = smc_em_bern_FoBaMo_v4(Sim,R,P)

%% makes code prettier :-)
B.sig2_c    = P.sigma_c^2*Sim.dt;
B.a         = Sim.dt/P.tau_c;
B.Aeff      = P.A*B.a;
B.kx        = P.k'*Sim.x;
if Sim.M>0
    B.sig2_h    = P.sigma_h.^2*Sim.dt;
    B.g         = 1-Sim.dt/P.tau_h;
    B.omega     = P.omega;
end
B.C_0       = P.C_0;
B.n         = P.n;
B.k_d       = P.k_d;
B.alpha     = P.alpha;
B.beta      = P.beta;
B.gamma     = P.gamma;
B.zeta      = P.zeta;

%% forward step
fprintf('\nTotal Number of Steps=%d\n',Sim.T)
if Sim.pf==0
    fprintf('\nprior forward step:          ')
    S = smc_em_bern_PriorSampl_v4(Sim,R,B);
elseif Sim.pf==1
    fprintf('\nk-mixture forward step:      ')
    S = smc_em_bern_k_comp_BackSampl_v3(Sim,R,B);
end

%% backward step
fprintf('\nbackward step:               ')
S.w_b = smc_em_bern_backwardPF_v9(Sim,S,B);
fprintf('\n')

%%   compute moments
M.nbar = sum(S.w_b.*S.n,1);
M.nvar = sum((repmat(M.nbar,Sim.N,1)-S.n).^2)/Sim.N;

M.Cbar = sum(S.w_b.*S.C,1);
M.Cvar = sum((repmat(M.Cbar,Sim.N,1)-S.C).^2)/Sim.N;

%%   compute mse
M.mse.n=sum(S.w_b.*(S.n-repmat(R.n,Sim.N,1)).^2);
M.mse.C=sum(S.w_b.*(S.C-repmat(R.C,Sim.N,1)).^2);
