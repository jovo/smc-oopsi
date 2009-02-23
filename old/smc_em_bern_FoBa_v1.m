function S = smc_em_bern_FoBa_v1(Sim,R,P)

%% makes code prettier :-)
B.sig2_c    = P.sigma_c^2*Sim.dt;
B.sig2_o    = P.sigma_o^2;
B.a         = 1-Sim.dt/P.tau_c;
B.beta      = P.beta;
B.kx        = P.k'*Sim.x;
if Sim.M>0
    B.sig2_h    = P.sigma_h.^2*Sim.dt;
    B.g         = 1-Sim.dt/P.tau_h;
    B.omega     = P.omega;
end

%% forward step
if Sim.van==false
    S = smc_em_bern_BackSampl_v1(Sim,R,B);
else
    S = smc_em_bern_PriorSampl_v1(Sim,R,B);
end

%% backward step
S.w_b = smc_em_bern_backwardPF_v6(Sim,S,B);