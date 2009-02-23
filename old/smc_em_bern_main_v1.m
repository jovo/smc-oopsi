function [B, S, M, E] = smc_em_bern_main_v1(Sim,R,P)

E(1:max(Sim.pfs),1:Sim.iters+1) = P;                                       %initialize parameter estimates

for pf=Sim.pfs
    for iter=1:Sim.iters

        %% makes code prettier :-)
        B.sig2_c    = E(pf,iter).sigma_c^2*Sim.dt;
        B.sig2_o    = E(pf,iter).sigma_o^2;
        B.a         = 1-Sim.dt/E(pf,iter).tau_c;
        B.beta      = E(pf,iter).beta;
        B.kx         = E(pf,iter).k'*Sim.x;
        if Sim.M>0
            B.sig2_h    = E(pf,iter).sigma_h.^2*Sim.dt;
            B.g         = 1-Sim.dt/E(pf,iter).tau_h;
            B.omega     = E(pf,iter).omega;
        end

        %% forward step
        Sim.van = false;
        if pf==1
            S{pf} = smc_em_bern_FohBack_v4c(Sim,R,B);
%         elseif pf==2
%             S{pf} = smc_em_bern_FohBack_v4d(Sim,R,B);
        else
            Sim.van = true;
            S{pf} = smc_em_bern_FoVanil_v1(Sim,R,B);
        end

        %% backward step
        S{pf}.w_b = smc_em_bern_backwardPF_v4(Sim,S{pf},B);

        %% get moments and mse's
        M(pf,iter) = smc_em_bern_moments_v1(Sim,S{pf},R);

        %% maximization step
        E(pf,iter+1) = E(pf,iter);
        %     E(pf+1iter+1) = smc_em_bern_mstep_v1(Sim,S,E(,pf+1iter));
        %     norm(E(pf+1,iter+1).k-P.k)
    end
end