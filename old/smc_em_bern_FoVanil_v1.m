function S = smc_em_bern_FoVanil_v1(Sim,R,B)

%% initialize stuff
S.p     = zeros(Sim.N,Sim.K);                   %initialize rate
S.I     = zeros(Sim.N,Sim.K);                   %initialize spike counts
S.C     = zeros(Sim.N,Sim.K);               %initialize calcium
S.h     = zeros(Sim.N,Sim.K,Sim.M);         %extize spike history terms
S.w_f   = 1/Sim.N*ones(Sim.N,Sim.K);          %initialize N_{eff}
S.w_b   = 1/Sim.N*ones(Sim.N,Sim.K);          %initialize N_{eff}
S.Neff  = Sim.N*ones(1,Sim.K_o);                %initialize N_{eff}

epsilon_c   = sqrt(B.sig2_c)*randn(Sim.N,Sim.K);%generate noise on c
U_sampl     = rand(Sim.N,Sim.K);                %random samples

if Sim.M>0  %if spike histories, generate noise on them
    S.h         = zeros(Sim.N,Sim.K,Sim.M);   %initialize spike history terms
    epsilon_h   = zeros(Sim.N, Sim.K, Sim.M);
    for m=1:Sim.M
        epsilon_h(:,:,m)   = sqrt(B.sig2_h(m))*randn(Sim.N,Sim.K);   %generate noise on h
    end
else	    %if no spike histories, generate p and sample I
    S.p(1,:) = 1-exp(-exp(B.kx)*Sim.dt);                                  %update rate for those particles with y_t<0
    S.p      = repmat(S.p(1,:),Sim.N,1);%make rate the same for each particle
    S.I      = U_sampl<S.p;             %generate random number to use for sampling
end

% preprocess stuff for stratified resampling
ints        = linspace(0,1,Sim.N+1);
diffs       = ints(2)-ints(1);
U_resamp    = repmat(ints(1:end-1),Sim.K_o,1)+diffs*rand(Sim.K_o,Sim.N);

%% loop-de-loop
for k=2:Sim.K

    % if h's, update h and I recursively
    if Sim.M>0                                      %update noise on h
        for m=1:Sim.M
            S.h(:,k,m)=B.g(m)*S.h(:,k-1,m)+S.I(:,k-1)+epsilon_h(:,k,m);
        end

        % update rate and sample spikes
        hs              = S.h(:,k,:);               %this is required for matlab to handle a m-by-n-by-p matrix
        h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
        y_t             = B.kx(k)+B.omega'*h';      %input to neuron
        S.p(:,k)        = 1-exp(-exp(y_t)*Sim.dt);  %update rate for those particles with y_t<0
        S.I(:,k)        = U_sampl(:,k)<S.p(:,k);    %sample
    end

    % sample C
    if mod(k,Sim.freq)==0 && B.sig2_o==0
        S.c(:,k)=R.C(k);
    else
        S.C(:,k)=B.a*S.C(:,k-1)+B.beta*S.I(:,k)+epsilon_c(:,k);
    end

    % stratified resample at every observation
    if mod(k,Sim.freq)==0
        if B.sig2_o==0
            S.w_f(:,k)  = 1/Sim.N*ones(Sim.N,1);
        else
            ln_w        = -0.5*(R.O(k)-S.C(:,k)).^2/B.sig2_o;     %compute log of weights
            ln_w        = ln_w-max(ln_w);                       %subtract the max to avoid rounding errors
            w           = exp(ln_w);                            %exponentiate to get actual weights
            S.w_f(:,k)  = w/sum(w);                             %normalize to define a legitimate distribution
        end        
        S = smc_em_bern_stratresamp_v5(Sim,S,k,U_resamp);
    end %resample
end