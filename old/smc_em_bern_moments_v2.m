function M = smc_em_bern_moments_v2(Sim,S,R)

%%   compute forward moments
M.fIbar = sum(S.w_f.*S.I,1);
M.fIvar = sum((repmat(M.fIbar,Sim.N,1)-S.I).^2)/Sim.N;

M.fCbar = sum(S.w_f.*S.C,1);
M.fCvar = sum((repmat(M.fCbar,Sim.N,1)-S.C).^2)/Sim.N;

%%   compute backward moments
M.bIbar = sum(S.w_b.*S.I,1);
M.bIvar = sum((repmat(M.bIbar,Sim.N,1)-S.I).^2)/Sim.N;

M.bCbar = sum(S.w_b.*S.C,1);
M.bCvar = sum((repmat(M.bCbar,Sim.N,1)-S.C).^2)/Sim.N;

%% compute MSE's
M.mse.I=sum(S.w_b.*(S.I-repmat(R.I,Sim.N,1)).^2);
M.mse.C=sum(S.w_b.*(S.C-repmat(R.C,Sim.N,1)).^2);