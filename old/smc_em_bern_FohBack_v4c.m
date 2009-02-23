function S = smc_em_bern_FohBack_v4c(Sim,R,B)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that Sim.M=0 (M is the # of spike history terms).
%
% The backward sampler has standard variance,  as approximated typically.
% Each particle still has its own backwards sampler,
% because spike histories make the p_k a function of the particle.
%
% this version only stores backwards sampling stuff as need.
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
S.C         = zeros(Sim.N,Sim.K);  %extize calcium
S.h         = zeros(Sim.N,Sim.K,Sim.M);         %extize spike history terms
S.w_f       = 1/Sim.N*ones(Sim.N,Sim.K);        %extize forward weights
S.w_b       = 1/Sim.N*ones(Sim.N,Sim.K);        %extize forward weights
S.Neff      = 1/Sim.N*ones(1,Sim.K_o);          %extize N_{eff}

% extize stuff needed for backwards sampling
O.mu_o      = zeros(Sim.N,Sim.freq);            %extize backwards mean
O.sig2_o    = zeros(Sim.N,Sim.freq);            %extize backwards variance

O.mu_I      = zeros(Sim.N,Sim.freq);            %extize mean spike indicator
O.sig2_I    = zeros(Sim.N,Sim.freq);            %extize variance spike indicator
O.sig2_c    = zeros(Sim.N,Sim.freq);            %extize variace of calcium

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
t               = Sim.freq;
O.mu_o(:,t)     = R.O(t);                       %initialize mean of P[O_s | C_s]
O.sig2_o(:,t)   = B.sig2_o;                     %initialize var of P[O_s | C_s]
O               = UpdateMoments(Sim,R,B,S,O,t); %recurse back to get P[O_s | C_s] before the first observation

%% loop-de-loop
for k=Sim.freq+1:Sim.K-Sim.freq

    %% sample h
    for m=1:Sim.M                               %for each spike history term
        S.h(:,k,m)=B.g(m)*S.h(:,k-1,m)+S.I(:,k-1)+epsilon_h(:,k,m);
    end

    %% sample I
    %  update p
    hs              = S.h(:,k,:);               %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_k             = B.kx(k)+B.omega'*h';      %input to neuron
    S.p(:,k)        = 1-exp(-exp(y_k)*Sim.dt);  %update p

    % compute P[I_k | h_k]
    ln_I    = [log(S.p(:,k)) log(1-S.p(:,k))];  %compute [log(spike) log(no spike)]

    % compute log G_n(I_k | O_s) for I_k=1 and I_k=0
    m1  = [B.a*S.C(:,k-1)+B.beta B.a*S.C(:,k-1)];%mean of P[C_k | C_{k-1}, I_k] for I_k=1 and I_k=0
    v1  = B.sig2_c;                             %var of P[C_k | C_{k-1}, I_k] for I_k=1 and I_k=0
    m2  = repmat(O.mu_o(:,k-t),1,2);            %mean of P[O_s | C_k] for I_k=1 and I_k=0
    v2  = repmat(O.sig2_o(:,k-t),1,2);          %var of P[O_s | C_k] for I_k=1 and I_k=0
    v   = (v1+v2);                              %var of G_n(I_k | O_s) for I_k=1 and I_k=0
    ln_G= -0.5*log(2*pi.*v)-.5*(m1-m2).^2./v;   %log G_n(I_k | O_s) for I_k=1 and I_k=0

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
    if B.sig2_o==0                              %if no observation noise
        S.C(:,k)=R.C(k);                        %at observations, C=O
        log_q_C = (1/Sim.N)*ones(Sim.N,1);      %likelihood of sampling that is negligible
    else                                        %if there is observation noise
        v_c     = (1./O.sig2_o(:,k-t)+1/B.sig2_c).^(-1);%variance of dist'n for sampling C
        m_c     = v_c.*(O.mu_o(:,k-t)./O.sig2_o(:,k-t)+(B.a*S.C(:,k-1)+B.beta*S.I(:,k))/B.sig2_c);%mean of dist'n for sampling C
        S.C(:,k)= normrnd(m_c,sqrt(v_c));       %sample C
        log_q_C = -0.5*(S.C(:,k)-m_c).^2./v_c;  %log prov of sampling the C_k that was sampled
    end
    %% update weights
    if mod(k,Sim.freq)==0            %at observations compute P(O|H)
        if B.sig2_o==0
            log_PO_H            = (1/Sim.N)*ones(Sim.N,1);
        else
            log_PO_H            = -0.5*(S.C(:,k)-R.O(k)).^2/B.sig2_o;
        end
    else                                        %when no observations are present, P[O|H^{(i)}] are all equal
        log_PO_H            = (1/Sim.N)*ones(Sim.N,1);
    end
    log_I           = ones(Sim.N,1);            %extize log sampling spikes
    log_I(sp)       = log(S.p(sp,k));           %compute log P(spike)
    log_I(nosp)     = log(1-S.p(nosp,k));       %compute log P(no spike)

    log_C_Cn        = -0.5*(S.C(:,k)-(B.a*S.C(:,k-1)+B.beta*S.I(:,k))).^2/B.sig2_c;%log P[C_k | C_{k-1}, I_k]

    log_q_n         = ones(Sim.N,1);            %initialize log q_n
    log_q_n(sp)     = log(q_n(sp,1));           %compute what was the log prob of sampling a spike
    log_q_n(nosp)   = log(1-q_n(nosp,1));       %or sampling no spike

    log_quotient    = log_PO_H + log_I + log_C_Cn - log_q_n - log_q_C;

    sum_logs        = log_quotient+log(S.w_f(:,k-1));   %update log(weights)
    w               = exp(sum_logs-max(sum_logs));      %exponentiate log(weights)
    S.w_f(:,k)      = w./sum(w);                        %normalize such that they sum to unity

    %% at observations
    if mod(k,Sim.freq)==0
%         figure(1), plot(R.C(1:k),'k','linewidth',2), hold on, plot(sum(S.w_f(:,1:k).*S.C(:,1:k)),'k'), plot(k,R.O(k),'ok')
        S = smc_em_bern_stratresamp_v5(Sim,S,k,U_resamp);   %stratified resample
        O = UpdateMoments(Sim,R,B,S,O,k);                   %estimate P[O_s | C_tt] for all k'<tt<s as a gaussian
        t = k;
    end
end %for time loop

end %function

%% update moments function
function O = UpdateMoments(Sim,R,B,S,O,k)

s               = Sim.freq;                     %find next observation time
O.mu_o(:,s)     = R.O(k+s);                     %initialize mean of P[O_s | C_s]
O.sig2_o(:,s)   = B.sig2_o;                     %initialize var of P[O_s | C_s]

hhat            = zeros(Sim.N,Sim.freq,Sim.M);  %extize hhat
phat            = zeros(Sim.N,Sim.freq+1);      %extize phat

hhat(:,1,:)     = S.h(:,k,:);                   %initialize hhat
phat(:,1)       = S.p(:,k);                     %initialize phat
for tt=1:s
    % update hhat
    for m=1:Sim.M                               %for each spike history term
        hhat(:,tt+1,m)=B.g(m)*hhat(:,tt,m)+phat(:,tt);
    end

    % update phat
    hs              = hhat(:,tt+1,:);           %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_t             = B.kx(tt+k)+B.omega'*h';   %input to neuron
    phat(:,tt+1)  = 1-exp(-exp(y_t)*Sim.dt);    %update phat
end

for tt=s:-1:2
    O.mu_I(:,tt)      = B.beta*phat(:,tt);                          %expected mean of spiking
    O.sig2_I(:,tt)    = B.beta^2*phat(:,tt).*(1-phat(:,tt));        %expected var of spiking

    O.sig2_c(:,tt)    = B.sig2_c+O.sig2_I(:,tt);                    %variance of C

    O.mu_o(:,tt-1)    = B.a^(-1)*(O.mu_o(:,tt)-O.mu_I(:,tt));       %mean of P[O_s | C_k]
    O.sig2_o(:,tt-1)  = B.a^(-2)*(O.sig2_c(:,tt)+O.sig2_o(:,tt));   %var of P[O_s | C_k]
end

end %function UpdateMoments