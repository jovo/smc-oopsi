function [S M] = smc_em_bern_main_v5(Sim,R,P)

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
    S = smc_em_bern_FohBack_v4c(Sim,R,B);
else
    S = smc_em_bern_FoVanil_v1(Sim,R,B);
end

%% backward step
S.w_b = smc_em_bern_backwardPF_v5(Sim,S,B);

%% get moments and mse's
M = smc_em_bern_moments_v2(Sim,S,R);

end %function