clear; clc;
load good_data2

Sim.van     = false;
P.sigma_o   = sqrt(P.sigma_c^2*Sim.dt);
sigma_o     = [0    10*P.sigma_o 20*P.sigma_o 50*P.sigma_o];
dtSample    = [     Sim.dt     5*Sim.dt    10*Sim.dt    20*Sim.dt];

for n=1:length(sigma_o)
    for i=1:length(dtSample)

        P.sigma_o   = sigma_o(n);
        epsilon_o   = P.sigma_o*randn(1,Sim.K);                   %generate noise on observations
        R.O         = R.C + epsilon_o;                        %add noise to observations
        Sim.dtSample= dtSample(i);
        Sim.K_o     = Sim.Nsec/Sim.dtSample;                %number of observations
        Sim.frac    = Sim.dtSample/Sim.dt;                  %# of samples per observation

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

        Sim.van= false;
        Sim.N  = 100;
        S{n,i} = smc_em_bern_FohBack_v4d(Sim,R,B);
%         Sim.van=true;
%         S{n,i} = smc_em_bern_FoVanil_v1(Sim,R,B);

        %% backward step
        S{n,i}.w_b = smc_em_bern_backwardPF_v4(Sim,S{n,i},B);

        %% get moments and mse's
        M(n,i) = smc_em_bern_moments_v1(Sim,S{n,i},R);

        Os{n,i} = R.O;
    end
end

GetArrayFig(Sim,R,M,Os,dtSample)