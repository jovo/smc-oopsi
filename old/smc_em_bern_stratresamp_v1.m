function S = smc_em_bern_stratresamp_v1(Sim,S,R,k,sig2_o,U_resamp)
% this function does stratified resampling for neuron models with spike
% history terms
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sim   : simulation parameters (eg, dt, N, etc.)
% S     : simulation states (eg, n, lambda, C, h)
% R     : "real" neuron's state
% k     : current time step index
% sig2_o: variance of observation noise
%
% i pre-compute a matrix of uniform (0,1) numbers using the following code
%     %preprocess stuff for stratified resampling
%     ints        = linspace(0,1,Sim.N+1);
%     diffs       = ints(2)-ints(1);
%     U_resamp    = repmat(ints(1:end-1),Sim.Nsamples,1)+diffs*rand(Sim.Nsamples,Sim.N);
%     Nresamp     = 0;
%
% and then i pass in U_resamp and Nresamp to this function.  this speeds
% things up so i don't need to generate random numbers with each sample
% time.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% S: this has particle states having resampled

% 1) re-sort to avoid bias due to ordering
new_ind=randperm(Sim.N);
S.p(:,k)=S.p(new_ind,k);
S.I(:,k)=S.I(new_ind,k);
S.C(:,k)=S.C(new_ind,k);
if Sim.M>0
    S.h(:,k,:)=S.h(new_ind,k,:);
end

% 2) compute cumulative weights
if sig2_o==0
    w = 1/Sim.N*ones(Sim.N,1);
else
    ln_w    = -0.5*(R.O(k)-S.C(:,k)).^2/sig2_o;     %compute log of weights
    ln_w    = ln_w-max(ln_w);                       %subtract the max to avoid rounding errors
    w       = exp(ln_w);                            %exponentiate to get actual weights
    w       = w/sum(w);                             %normalize to define a legitimate distribution
end
cum_w   = cumsum(w);                            %get cumulative sum of weights for sampling purposese
S.Neff(k/Sim.frac)  = 1/sum(w.^2);              %store N_{eff}

% 3) resample
Nresamp=k/Sim.frac;                             %increase sample counter
for n=1:Sim.N                                   %for each particle
    Nparticle   = find(U_resamp(Nresamp,n)<cum_w,1);%sample
    S.p(n,k)    = S.p(Nparticle,k);             %reset rate
    S.I(n,k)    = S.I(Nparticle,k);             %reset n
    S.C(n,k)    = S.C(Nparticle,k);             %reset C
    if Sim.M>0                                  %if spike history terms
        S.h(n,k,:) = S.h(Nparticle,k,:);        %reset all h's
    end
end