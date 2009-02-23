function S = strat_resample(Sim,S,t,U_resamp)
% this function does stratified resampling for our model
%
% Inputs:
% Sim       : simulation parameters (eg, dt, N, etc.)
% S         : simulation states (eg, n, lambda, C, h)
% t         : current time step index
% U_resamp  : resampling matrix (so that i needed generate random numbers
%             which each resample, but rather, can just call them

% Outputs:
% S: this has particle states having resampled them, and equalized weights

Nresamp=t/Sim.freq;                         % increase sample counter
S.Neff(Nresamp)  = 1/sum(S.w_f(:,t).^2);    % store N_{eff}

% if weights are degenerate or we are doing prior sampling then resample
if S.Neff(Nresamp) < Sim.N/2 || Sim.pf==0
    [foo,ind]           = histc(U_resamp(Nresamp,:),[0  cumsum(S.w_f(:,t))']);
    [ri,ri]             = sort(rand(Sim.N,1));      % these 3 lines stratified resample
    ind                 = ind(ri);
    S.p(:,t-Sim.freq:t) = S.p(ind,t-Sim.freq:t);    % resample probabilities (necessary?)
    S.n(:,t-Sim.freq:t) = S.n(ind,t-Sim.freq:t);    % resample calcium
    S.C(:,t-Sim.freq:t) = S.C(ind,t-Sim.freq:t);    % resample calcium
    if Sim.M>0                                      % if spike history terms
        S.h(:,t-Sim.freq:t,:) = S.h(ind,t-Sim.freq:t,:);% resample all h's
    end
    S.w_f(:,t-Sim.freq:t)= 1/Sim.N*ones(Sim.N,1);   % reset weights
end %function