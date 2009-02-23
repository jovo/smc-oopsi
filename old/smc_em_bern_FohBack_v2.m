function S = smc_em_bern_FohBack_v2(Sim,R,B)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that Sim.M=0 (M is the # of spike history terms).  The
% backward sampler has increased variance by making the var of P[O_s | C_k] 
% an order of magnitude bigger. each particle still has its own backwards
% sampler though, because spike histories make the p_k a function
% of the particle.
% 
%% inputs: 
% Sim:  simulation parameters
% R:    "real" neuron data
% B:    current parameter estimates
% 
%% outputs: 
% S:    simulation states

%% initialize stuff
% initialize particle info
S.p         = zeros(Sim.N,Sim.K);               %initialize rate
S.I         = zeros(Sim.N,Sim.K);               %initialize spike counts
S.C         = R.O(Sim.frac)*ones(Sim.N,Sim.K);  %initialize calcium
S.h         = 0*ones(Sim.N,Sim.K,Sim.M);   %initialize spike history terms
S.w_f       = 1/Sim.N*ones(Sim.N,Sim.K);        %initialize forward weights
S.Neff      = 1/Sim.N*ones(1,Sim.K_o);          %initialize N_{eff}

% initialize stuff needed for backwards sampling
S.mu_o      = zeros(Sim.N,Sim.K);           %initialize backwards mean
S.sig2_o    = zeros(Sim.N,Sim.K);           %initialize backwards variance

S.mu_I      = zeros(Sim.N,Sim.K);           %initialize mean spike indicator
S.sig2_I    = zeros(Sim.N,Sim.K);           %initialize variance spike indicator
S.sig2_c    = zeros(Sim.N,Sim.K);           %initialize variace of calcium

% preprocess stuff for stratified resampling
ints        = linspace(0,1,Sim.N+1);    %generate intervals
diffs       = ints(2)-ints(1);          %generate interval size
U_resamp    = repmat(ints(1:end-1),Sim.K_o,1)+diffs*rand(Sim.K_o,Sim.N); %resampling matrix

% initialize misc stuff
U_sampl     = rand(Sim.N,Sim.K);        %generate random number to use for sampling

% generate noise on h
epsilon_h   = zeros(Sim.N, Sim.K, Sim.M);       %initialize noise on h
for m=1:Sim.M                                   %generate noise on h
    epsilon_h(:,:,m)   = sqrt(B.sig2_h(m))*randn(Sim.N,Sim.K);
end

%% take this out, just for debugging
S.p2(1,:)    = 1-exp(-exp(B.kx));                  %update rate for those particles with s<0

% initialize backwards distributions
S.mu_o(:,Sim.frac)  = R.O(Sim.frac);                %expected mean of O
S.sig2_o(:,Sim.frac)= B.sig2_o;                     %expected var of O
S                   = UpdateMomens(Sim,R,B,S,Sim.frac);
    
%% loop-de-loop
for k=Sim.frac:Sim.K-Sim.frac

    %% sample h
    for m=1:Sim.M
        S.h(:,k,m)=B.g(m)*S.h(:,k-1,m)+S.I(:,k-1)+epsilon_h(:,k,m);
    end

    %% update p
    hs              = S.h(:,k,:);               %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_t             = B.kx(k)+B.omega'*h';        %input to neuron
    S.p(:,k)        = 1-exp(-exp(y_t)*Sim.dt);         %update p

    %% sample I
    ln_I    = [log(S.p(:,k)) log(1-S.p(:,k))];          %compute [log(spike) log(no spike)]
    mu_ck   = [B.a*S.C(:,k-1)+B.beta B.a*S.C(:,k-1)];   %mean when spiked and not spiked

    sig2_zk = repmat((1/B.sig2_c + 1./S.sig2_o(:,k)).^(-1),1,2);        %update sig2_zk
    mu_zk   = sig2_zk.*(mu_ck/B.sig2_c + repmat(S.mu_o(:,k)./S.sig2_o(:,k),1,2));                     %update mu_zk
    ln_Z    = log(sig2_zk./sqrt(2*pi*repmat(S.sig2_o(:,k),1,2)*B.sig2_c)) + 0.5*(mu_zk.^2./sig2_zk-repmat(S.mu_o(:,k).^2./S.sig2_o(:,k),1,2)-mu_ck.^2/B.sig2_c);  %update ln_z (ignore constant)
    ln_q_n  = ln_I + ln_Z;                              %log of sampling dist
    mx      = max(ln_q_n,[],2);                         %find max of each column
    mx2     = repmat(mx,1,2);                           %matricize
    q_n     = exp(ln_q_n-mx2);                          %subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
    q_n     = q_n./repmat(sum(q_n,2),1,2);              %normalize to make this a true sampling distribution (ie, sums to 1)
    S.I(:,k)= U_sampl(:,k)<q_n(:,1);                    %sample I
    sp      = S.I(:,k)==1;
    nosp    = S.I(:,k)==0;

    %% sample C
    v_c     = (1./S.sig2_o(:,k)+1/B.sig2_c).^(-1);          %variance of dist'n for sampling C
    m_c     = v_c.*(S.mu_o(:,k)./S.sig2_o(:,k)+(B.a*S.C(:,k-1)+B.beta*S.I(:,k))/B.sig2_c);   %mean of dist'n for sampling C
    S.C(:,k)= normrnd(m_c,sqrt(v_c));                   %sample C

    %% update weights
    if mod(k*Sim.dt,Sim.dtSample)==0                    %at observations compute P(O|H)
        log_PO_H            = -0.5*(S.C(:,k)-R.O(k)).^2/B.sig2_o;
    else
        log_PO_H            = (1/Sim.N)*ones(Sim.N,1);
    end
    log_I           = ones(Sim.N,1);
    log_I(sp)       = log(S.p(sp,k));          %compute [log(spike) log(no spike)]
    log_I(nosp)     = log(1-S.p(nosp,k));          %compute [log(spike) log(no spike)]

    log_C_Cn        = -0.5*(S.C(:,k)-(B.a*S.C(:,k-1)+B.beta*S.I(:,k))).^2/B.sig2_c;

    log_q_n         = ones(Sim.N,1);
    log_q_n(sp)     = log(q_n(sp,1));           %compute [log(spike) log(no spike)]
    log_q_n(nosp)   = log(1-q_n(nosp,1));         %compute [log(spike) log(no spike)]

    log_q_C         = -0.5*(S.C(:,k)-m_c).^2./v_c;

    log_quotient    = log_PO_H + log_I + log_C_Cn - log_q_n - log_q_C;

    sum_logs        = log_quotient+log(S.w_f(:,k-1));                    %update log(weights)
    w               = exp(sum_logs-max(sum_logs));
    S.w_f(:,k)      = w./sum(w);                                       %normalize such that they sum to unity

    %% at observations
    if mod(k*Sim.dt,Sim.dtSample)==0
        S = smc_em_bern_stratresamp_v2(Sim,S,k,U_resamp);   %stratified resample
        S = UpdateMomens(Sim,R,B,S,k);                      %estimate P[O_s | C_tt] for all k'<tt<s as a gaussian
    end
end %for time loop

end %function

%% update moments function
function S = UpdateMomens(Sim,R,B,S,k)

s               = k+Sim.frac;
S.mu_o(:,s)     = R.O(s);
S.sig2_o(:,s)   = B.sig2_o;

hhat            = zeros(Sim.N,Sim.frac,Sim.M);
phat            = zeros(Sim.N,Sim.frac+1);

hhat(:,1,:)     = S.h(:,k,:);
phat(:,1)       = S.p(:,k);
phat2           = S.p(:,k:s);
for tt=k+1:s
    %% update hhat
    for m=1:Sim.M
        hhat(:,tt-k+1,m)=B.g(m)*hhat(:,tt-k,m)+phat(:,tt-k);
    end
    %% update phat
    hs              = hhat(:,tt-k+1,:);         %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_t             = B.kx(tt)+B.omega'*h';      %input to neuron
    phat(:,tt-k+1)  = 1-exp(-exp(y_t)*Sim.dt);         %update p
end

for tt=s:-1:k+2
    S.mu_I(:,tt)      = B.beta*phat(:,tt-k);                     %expected mean of spiking
    S.sig2_I(:,tt)    = B.beta^2*phat(:,tt-k).*(1-phat(:,tt-k));     %expected var of spiking

    S.sig2_c(:,tt)    = B.sig2_c+S.sig2_I(:,tt);                %variance of C

    S.mu_o(:,tt-1)    = B.a^(-1)*(S.mu_o(:,tt)-S.mu_I(:,tt));     %mean of P[O_s | C_k]
    S.sig2_o(:,tt-1)  = 10*B.a^(-2)*(S.sig2_c(:,tt)+S.sig2_o(:,tt)); %var of P[O_s | C_k]
end
S.mu_I(:,k+1)      = B.beta*phat(:,1);
S.sig2_I(:,k+1)    = B.beta^2*phat(:,1).*(1-phat(:,1));

S.sig2_c(:,k+1)    = B.sig2_c+S.sig2_I(:,k+1);

end %function UpdateMomens