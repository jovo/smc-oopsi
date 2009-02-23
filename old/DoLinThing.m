function [S M] = DoLinThing(Sim,R)

P.k         = 1;                %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = .1;               %std of noise on calcium
P.sigma_o   = 10*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = -1;                %jump size for h after spike
P.tau_h     = 1;                %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

% do EM recursion
Sim.pf  = 3;            %not vanilla particle filtering
R.O     = (R.F-1)/2;
[S M]   = smc_em_bern_FoBaMo_v3(Sim,R,P);%do forward-backward and get moments
