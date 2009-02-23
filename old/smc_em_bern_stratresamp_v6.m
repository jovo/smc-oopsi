function S = smc_em_bern_stratresamp_v6(Sim,S,t,U_resamp)
% this function does stratified resampling for neuron models with spike
% history terms
%
%% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sim       : simulation parameters (eg, dt, N, etc.)
% S         : simulation states (eg, n, lambda, C, h)
% t         : current time step index
% U_resamp  : resampling matrix (so that i needed generate random numbers
%             which each resample, but rather, can just call them

%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% S: this has particle states having resampled them, and equalized weights

%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nresamp=t/Sim.freq;                             %increase sample counter
S.Neff(Nresamp)  = 1/sum(S.w_f(:,t).^2);              %store N_{eff}

if S.Neff(Nresamp) < Sim.N/2

    [foo,ind]   = histc(U_resamp(Nresamp,:),[0  cumsum(S.w_f(:,t))']);
    [ri,ri]     = sort(rand(Sim.N,1));
    ind         = ind(ri);
    S.p(:,t)    = S.p(ind,t);
    S.n(:,t)    = S.n(ind,t);
    S.C(:,t)    = S.C(ind,t);
    if Sim.M>0                                  %if spike history terms
        S.h(:,t,:) = S.h(ind,t,:);        %resample all h's
    end
    S.w_f(:,t)=1/Sim.N*ones(Sim.N,1);             %reset weights

end %if Neff<N/2

end %function