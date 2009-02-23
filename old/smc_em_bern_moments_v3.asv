function M = smc_em_bern_moments_v3(Sim,S)

%%   compute moments
M.nbar = sum(S.w_b.*S.n,1);
M.nvar = sum((repmat(M.nbar,Sim.N,1)-S.n).^2)/Sim.N;

M.Cbar = sum(S.w_b.*S.C,1);
M.Cvar = sum((repmat(M.Cbar,Sim.N,1)-S.C).^2)/Sim.N;