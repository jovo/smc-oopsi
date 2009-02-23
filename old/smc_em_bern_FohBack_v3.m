function S = smc_em_bern_FohBack_v3(Sim,R,B)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that Sim.M=0 (M is the # of spike history terms).  
% 
% The backward sampler has standard variance,  as approximated typically.
% Each particle still has its own backwards sampler, 
% because spike histories make the p_k a function of the particle.
% 
% inputs: 
% Sim:  simulation parameters
% R:    "real" neuron data
% B:    current parameter estimates
% 
% outputs: 
% S:    simulation states

%% allocate memory and initialize stuff
% extize particle info
S.p         = zeros(Sim.N,Sim.K);               %extize rate
S.I         = zeros(Sim.N,Sim.K);               %extize spike counts
S.C         = R.O(Sim.frac)*ones(Sim.N,Sim.K);  %extize calcium
S.h         = zeros(Sim.N,Sim.K,Sim.M);         %extize spike history terms
S.w_f       = 1/Sim.N*ones(Sim.N,Sim.K);        %extize forward weights
S.Neff      = 1/Sim.N*ones(1,Sim.K_o);          %extize N_{eff}

% extize stuff needed for backwards sampling
S.mu_o      = zeros(Sim.N,Sim.K);               %extize backwards mean
S.sig2_o    = zeros(Sim.N,Sim.K);               %extize backwards variance

S.mu_I      = zeros(Sim.N,Sim.K);               %extize mean spike indicator
S.sig2_I    = zeros(Sim.N,Sim.K);               %extize variance spike indicator
S.sig2_c    = zeros(Sim.N,Sim.K);               %extize variace of calcium

% preprocess stuff for stratified resampling
ints        = linspace(0,1,Sim.N+1);            %generate intervals
diffs       = ints(2)-ints(1);                  %generate interval size
U_resamp    = repmat(ints(1:end-1),Sim.K_o,1)+diffs*rand(Sim.K_o,Sim.N); %resampling matrix

% extize misc stuff
U_sampl     = rand(Sim.N,Sim.K);                %generate random number to use for sampling

% generate noise on h
epsilon_h   = zeros(Sim.N, Sim.K, Sim.M);       %generate noise on h
for m=1:Sim.M                                   %add noise to each h
    epsilon_h(:,:,m)   = sqrt(B.sig2_h(m))*randn(Sim.N,Sim.K);
end

% initialize backwards distributions
S.mu_o(:,Sim.frac)  = R.O(Sim.frac);            %initialize mean of P[O_s | C_s]
S.sig2_o(:,Sim.frac)= B.sig2_o;                 %initialize var of P[O_s | C_s]
S                   = UpdateMomens(Sim,R,B,S,Sim.frac);%recurse back to get P[O_s | C_s] before the first observation
    
%% loop-de-loop
for k=Sim.frac:Sim.K-Sim.frac

    %% sample h
    for m=1:Sim.M                               %for each spike history term
        S.h(:,k,m)=B.g(m)*S.h(:,k-1,m)+S.I(:,k-1)+epsilon_h(:,k,m);
    end

    %% sample I
    %  update p
    hs              = S.h(:,k,:);               %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_t             = B.kx(k)+B.omega'*h';      %input to neuron
    S.p(:,k)        = 1-exp(-exp(y_t)*Sim.dt);  %update p

    % compute P[I_k | h_k]
    ln_I    = [log(S.p(:,k)) log(1-S.p(:,k))];  %compute [log(spike) log(no spike)]

    % compute log G_n(I_k | O_s) for I_k=1 and I_k=0
    m1  = [B.a*S.C(:,k-1)+B.beta B.a*S.C(:,k-1)];%mean of P[C_k | C_{k-1}, I_k] for I_k=1 and I_k=0 
    v1  = B.sig2_c;                             %var of P[C_k | C_{k-1}, I_k] for I_k=1 and I_k=0
    m2  = repmat(S.mu_o(:,k),1,2);              %mean of P[O_s | C_k] for I_k=1 and I_k=0
    v2  = repmat(S.sig2_o(:,k),1,2);            %var of P[O_s | C_k] for I_k=1 and I_k=0
    v   = v1+v2;                                %var of G_n(I_k | O_s) for I_k=1 and I_k=0
    ln_G= -0.5*log(2*pi.*v)-.5*(m1-m2).^2./v; %log G_n(I_k | O_s) for I_k=1 and I_k=0

    % compute q(I_k | h_k, O_s)
    ln_q_n  = ln_I + ln_G;                      %log of sampling dist
    mx      = max(ln_q_n,[],2);                 %find max of each column
    mx2     = repmat(mx,1,2);                   %matricize
    q_n     = exp(ln_q_n-mx2);                  %subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
    q_n     = q_n./repmat(sum(q_n,2),1,2);      %normalize to make this a true sampling distribution (ie, sums to 1)
        
    % sample I
    S.I(:,k)= U_sampl(:,k)<q_n(:,1);            %sample I        
    sp      = S.I(:,k)==1;                      %store index of which samples spiked
    nosp    = S.I(:,k)==0;                      %and which did not

    %% sample C
    v_c     = (1./S.sig2_o(:,k)+1/B.sig2_c).^(-1);%variance of dist'n for sampling C
    m_c     = v_c.*(S.mu_o(:,k)./S.sig2_o(:,k)+(B.a*S.C(:,k-1)+B.beta*S.I(:,k))/B.sig2_c);%mean of dist'n for sampling C
    S.C(:,k)= normrnd(m_c,sqrt(v_c));           %sample C

    %% update weights
    if mod(k*Sim.dt,Sim.dtSample)==0            %at observations compute P(O|H)
        log_PO_H            = -0.5*(S.C(:,k)-R.O(k)).^2/B.sig2_o;
    else                                        %when no observations are present, P[O|H^{(i)}] are all equal
        log_PO_H            = (1/Sim.N)*ones(Sim.N,1);
    end
    log_I           = ones(Sim.N,1);
    log_I(sp)       = log(S.p(sp,k));           %compute log P(spike)
    log_I(nosp)     = log(1-S.p(nosp,k));       %compute log P(no spike)

    log_C_Cn        = -0.5*(S.C(:,k)-(B.a*S.C(:,k-1)+B.beta*S.I(:,k))).^2/B.sig2_c;%log P[C_k | C_{k-1}, I_k]

    log_q_n         = ones(Sim.N,1);            %initialize log q_n
    log_q_n(sp)     = log(q_n(sp,1));           %compute what was the lo prob of sampling a spike
    log_q_n(nosp)   = log(1-q_n(nosp,1));       %or sampling no spike

    log_q_C         = -0.5*(S.C(:,k)-m_c).^2./v_c;%log prov of sampling the C_k that was sampled

    log_quotient    = log_PO_H + log_I + log_C_Cn - log_q_n - log_q_C;

    sum_logs        = log_quotient+log(S.w_f(:,k-1));%update log(weights)
    w               = exp(sum_logs-max(sum_logs));%exponentiate log(weights)
    S.w_f(:,k)      = w./sum(w);                %normalize such that they sum to unity

    %% at observations
    if mod(k*Sim.dt,Sim.dtSample)==0
        S = smc_em_bern_stratresamp_v2(Sim,S,k,U_resamp);   %stratified resample
        S = UpdateMomens(Sim,R,B,S,k);                      %estimate P[O_s | C_tt] for all k'<tt<s as a gaussian
    end
end %for time loop

end %function

%% update moments function
function S = UpdateMomens(Sim,R,B,S,k)

s               = k+Sim.frac;                   %find next observation time
S.mu_o(:,s)     = R.O(s);                       %initialize mean of P[O_s | C_s]
S.sig2_o(:,s)   = B.sig2_o;                     %initialize var of P[O_s | C_s]

hhat            = zeros(Sim.N,Sim.frac,Sim.M);  %extize hhat
phat            = zeros(Sim.N,Sim.frac+1);      %extize phat

hhat(:,1,:)     = S.h(:,k,:);                   %initialize hhat
phat(:,1)       = S.p(:,k);                     %initialize phat
for tt=k+1:s
    % update hhat
    for m=1:Sim.M                               %for each spike history term
        hhat(:,tt-k+1,m)=B.g(m)*hhat(:,tt-k,m)+phat(:,tt-k);
    end
    
    % update phat
    hs              = hhat(:,tt-k+1,:);         %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_t             = B.kx(tt)+B.omega'*h';     %input to neuron
    phat(:,tt-k+1)  = 1-exp(-exp(y_t)*Sim.dt);  %update phat
end

for tt=s:-1:k+2
    S.mu_I(:,tt)      = B.beta*phat(:,tt-k);                        %expected mean of spiking
    S.sig2_I(:,tt)    = B.beta^2*phat(:,tt-k).*(1-phat(:,tt-k));    %expected var of spiking

    S.sig2_c(:,tt)    = B.sig2_c+S.sig2_I(:,tt);                    %variance of C

    S.mu_o(:,tt-1)    = B.a^(-1)*(S.mu_o(:,tt)-S.mu_I(:,tt));       %mean of P[O_s | C_k]
    S.sig2_o(:,tt-1)  = B.a^(-2)*(S.sig2_c(:,tt)+S.sig2_o(:,tt));   %var of P[O_s | C_k]
end

end %function UpdateMomens