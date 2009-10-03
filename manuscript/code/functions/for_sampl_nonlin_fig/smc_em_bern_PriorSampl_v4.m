function S = smc_em_bern_PriorSampl_v4(Sim,R,B)
% this code implements the prior sampler for our model
% 
% Input:
% Sim:  simulation parameters
% R:    real neuron data
% B:    parameter estimates
% 
% Output:
% S:    values of particles


%% initialize stuff
S.p     = zeros(Sim.N,Sim.T);                   %initialize rate
S.n     = zeros(Sim.N,Sim.T);                   %initialize spike counts
S.C     = B.C_init*ones(Sim.N,Sim.T);           %initialize calcium
S.w_f   = 1/Sim.N*ones(Sim.N,Sim.T);            %initialize forward weights
S.w_b   = 1/Sim.N*ones(Sim.N,Sim.T);            %initialize backward weights
S.Neff  = Sim.N*ones(1,Sim.T_o);                %initialize N_{eff}

A.epsilon_c     = sqrt(B.sig2_c)*randn(Sim.N,Sim.T);%generate noise on C
A.U_sampl       = rand(Sim.N,Sim.T);            %generate matrix to sample from

if Sim.M>0                                      %if spike histories, generate noise on them
    S.h         = zeros(Sim.N,Sim.T,Sim.M);     %initialize spike history terms
    A.epsilon_h = zeros(Sim.N, Sim.T, Sim.M);   %generate noise on h
    for m=1:Sim.M                               %for each h_m: m=1,...,K
        A.epsilon_h(:,:,m) = sqrt(B.sig2_h(m))*randn(Sim.N,Sim.T);
    end
else                                            %if no spike histories, generate p and sample I
    S.p(1,:) = 1-exp(-exp(B.kx)*Sim.dt);        %update rate 1 particle
    S.p      = repmat(S.p(1,:),Sim.N,1);        %make rate the same for each particle (as they are independent of spike history, p(spike) is identical for each
    S.n      = A.U_sampl<S.p;                   %sample spikes for all particles and all time
end

% preprocess stuff for stratified resampling
ints        = linspace(0,1,Sim.N+1);            %make a vector [0:1/N, 1/N:2/N,..., (N-1)/N:1]
diffs       = ints(2)-ints(1);                  %get the size of each increment, ie, 1/N
A.U_resamp  = repmat(ints(1:end-1),Sim.T_o,1)+diffs*rand(Sim.T_o,Sim.N);%generate a matrix from which we stratified resample

%% loop-de-loop
for t=2:Sim.T

    if Sim.M>0                                      %update noise on h
        for m=1:Sim.M                               %for each h_m
            S.h(:,t,m)=B.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
        end

        % update rate and sample spikes
        hs              = S.h(:,t,:);               %this is required for matlab to handle a m-by-n-by-p matrix
        h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
        y_t             = B.kx(t)+B.omega'*h';      %input to neuron
        S.p(:,t)        = 1-exp(-exp(y_t)*Sim.dt);  %update rate for particles
    end
    S.n(:,t)            = A.U_sampl(:,t)<S.p(:,t);  %sample
    S.C(:,t)        = (1-B.a)*S.C(:,t-1)+B.A*S.n(:,t)+B.a*B.C_0+A.epsilon_c(:,t);%sample C

    % stratified resample at every observation
    if mod(t,Sim.freq)==0
        F_mu        = B.alpha*Hill_v1(B,S.C(:,t))+B.beta;   %compute E[F_t]
        F_var       = B.gamma*Hill_v1(B,S.C(:,t)) + B.zeta; %compute V[F_t]
        ln_w        = -0.5*(R.F(t)-F_mu).^2./F_var + log(F_var)/2;%compute log of weights
        ln_w        = ln_w-max(ln_w);                       %subtract the max to avoid rounding errors
        w           = exp(ln_w);                            %exponentiate to get actual weights
        S.w_f(:,t)  = w/sum(w');                            %normalize to define a legitimate distribution
        S = smc_em_bern_stratresamp_v9(Sim,S,t,A.U_resamp);
    end

    if mod(t,100)==0
        if t<1000              %print # of observations
            fprintf('\b\b\b%d',t)
        elseif t<10000
            fprintf('\b\b\b\b%d',t)
        elseif t<100000
            fprintf('\b\b\b\b\b%d',t)
        end
    end
end %lood-de-loop

end %prior sampler