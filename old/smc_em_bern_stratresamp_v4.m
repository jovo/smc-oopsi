function S = smc_em_bern_stratresamp_v4(Sim,S,k,U_resamp)
% this function does stratified resampling for neuron models with spike
% history terms
%
%% inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Sim       : simulation parameters (eg, dt, N, etc.)
% S         : simulation states (eg, n, lambda, C, h)
% k         : current time step index
% U_resamp  : resampling matrix (so that i needed generate random numbers
%             which each resample, but rather, can just call them

%% outputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% S: this has particle states having resampled them, and equalized weights

%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nresamp=k/Sim.frac;                             %increase sample counter
S.Neff(Nresamp)  = 1/sum(S.w_f(:,k).^2);        %Compute N_{eff}

if S.Neff(Nresamp) < Sim.N/2

    [foo,ind]   = histc(U_resamp(Nresamp,:),[0  cumsum(S.w_f(:,k))']);
    [ri,ri]     = sort(rand(Sim.N,1));
    ind         = ind(ri);
    S.p(:,k)    = S.p(ind,k);
    S.I(:,k)    = S.I(ind,k);
    S.C(:,k)    = S.C(ind,k);
    if Sim.M>0                                  %if spike history terms
        S.h(:,k,:) = S.h(ind,k,:);        %resample all h's
    end
    S.w_f(:,k)=1/Sim.N*ones(Sim.N,1);             %reset weights

end %if Neff<N/2

end %function