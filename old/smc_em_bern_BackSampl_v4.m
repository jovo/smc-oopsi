function S = smc_em_bern_BackSampl_v4(Sim,R,B)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that Sim.M=0 (M is the # of spike history terms).
%
% The backward sampler has standard variance,  as approximated typically.
% Each particle has the SAME backwards sampler, initialized at E[h_t]
% this only does spike history stuff if necessary
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
S.n         = zeros(Sim.N,Sim.T);               %extize spike counts
S.C         = zeros(Sim.N,Sim.T);               %extize calcium
S.w_f       = 1/Sim.N*ones(Sim.N,Sim.T);        %extize forward weights
S.w_b       = 1/Sim.N*ones(Sim.N,Sim.T);        %extize forward weights
S.Neff      = 1/Sim.N*ones(1,Sim.T_o);          %extize N_{eff}

% extize stuff needed for backwards sampling
O.mu_o      = zeros(1,Sim.freq);                %extize backwards mean
O.sig2_o    = zeros(1,Sim.freq);                %extize backwards variance

O.mu_n      = zeros(1,Sim.freq);                %extize mean spike indicator
O.sig2_n    = zeros(1,Sim.freq);                %extize variance spike indicator
O.sig2_c    = zeros(1,Sim.freq);                %extize variace of calcium

% preprocess stuff for stratified resampling
ints        = linspace(0,1,Sim.N+1);            %generate intervals
diffs       = ints(2)-ints(1);                  %generate interval size
U_resamp    = repmat(ints(1:end-1),Sim.T_o,1)+diffs*rand(Sim.T_o,Sim.N); %resampling matrix

% extize misc stuff
U_sampl     = rand(Sim.N,Sim.T);                %generate random number to use for sampling
oney        = ones(Sim.N,1);                    %a vector of ones to call for various purposes to speed things up

% if spike histories
if Sim.M>0
    S.p         = zeros(Sim.N,Sim.T);           %extize prob(spike)
    S.h         = zeros(Sim.N,Sim.T,Sim.M);     %extize spike history terms
    epsilon_h   = zeros(Sim.N, Sim.T, Sim.M);   %generate noise on h
    for m=1:Sim.M                               %add noise to each h
        epsilon_h(:,:,m)   = sqrt(B.sig2_h(m))*randn(Sim.N,Sim.T);
    end
    % if not, comput P[n_t] for all t
else
    S.p         = 1-exp(-exp(B.kx)*Sim.dt);
end
% initialize backwards distributions
ss              = Sim.freq;
O.mu_o(:,ss)    = R.O(ss);                      %initialize mean of P[O_s | C_s]
O.sig2_o(:,ss)  = B.sig2_o;                     %initialize var of P[O_s | C_s]
O               = UpdateMoments(Sim,R,B,S,O,ss);%recurse back to get P[O_s | C_s] before the first observation

%% loop-de-loop
for t=Sim.freq+1:Sim.T-Sim.freq

    % if spike histories, sample h and update p
    if Sim.M>0                                  %update noise on h
        %% sample h
        for m=1:Sim.M                           %for each spike history term
            S.h(:,t,m)=B.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+epsilon_h(:,t,m);
        end

        %  update p
        hs              = S.h(:,t,:);           %this is required for matlab to handle a m-by-n-by-p matrix
        h(:,1:Sim.M)    = hs(:,1,1:Sim.M);      %this too
        S.p(:,t)        = 1-exp(-exp(B.kx(t)+B.omega'*h')*Sim.dt);
    end

    % compute P[n_k | h_k]
    ln_n    = [log(S.p(:,t)) log(1-S.p(:,t))];  %compute [log(spike) log(no spike)]

    % compute log G_n(n_k | O_s) for n_k=1 and n_k=0
    m1  = [B.a*S.C(:,t-1)+B.beta B.a*S.C(:,t-1)];%mean of P[C_k | C_{t-1}, n_k] for n_k=1 and n_k=0
    v1  = B.sig2_c;                             %var of P[C_k | C_{t-1}, n_k] for n_k=1 and n_k=0
    m2  = repmat(O.mu_o(t-ss),Sim.N,2);         %mean of P[O_s | C_k] for n_k=1 and n_k=0
    v2  = repmat(O.sig2_o(t-ss),Sim.N,2);       %var of P[O_s | C_k] for n_k=1 and n_k=0
    v   = (v1+v2);                              %var of G_n(n_k | O_s) for n_k=1 and n_k=0
    ln_G= -0.5*log(2*pi.*v)-.5*(m1-m2).^2./v;   %log G_n(n_k | O_s) for n_k=1 and n_k=0

    % compute q(n_k | h_k, O_s)
    ln_q_n  = ln_n + ln_G;                      %log of sampling dist
    mx      = max(ln_q_n,[],2);                 %find max of each column
    mx2     = repmat(mx,1,2);                   %matricize
    q_n     = exp(ln_q_n-mx2);                  %subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
    q_n     = q_n./repmat(sum(q_n,2),1,2);      %normalize to make this a true sampling distribution (ie, sums to 1)

    % sample n
    S.n(:,t)= U_sampl(:,t)<q_n(:,1);            %sample n
    sp      = S.n(:,t)==1;                      %store index of which samples spiked
    nosp    = S.n(:,t)==0;                      %and which did not

    %% sample C
    if B.sig2_o==0                              %if no observation noise
        S.C(:,t)=R.C(t);                        %at observations, C=O
        log_q_C = (1/Sim.N)*oney;               %likelihood of sampling that is negligible
    else                                        %if there is observation noise
        v_c     = (1./O.sig2_o(t-ss)+1/B.sig2_c).^(-1);%variance of dist'n for sampling C
        m_c     = v_c.*(O.mu_o(t-ss)./O.sig2_o(t-ss)+(B.a*S.C(:,t-1)+B.beta*S.n(:,t))/B.sig2_c);%mean of dist'n for sampling C
        S.C(:,t)= normrnd(m_c,sqrt(v_c));       %sample C
        log_q_C = -0.5*(S.C(:,t)-m_c).^2./v_c;  %log prov of sampling the C_k that was sampled
    end
    %% update weights
    if mod(t,Sim.freq)==0                       %at observations compute P(O|H)
        if B.sig2_o==0                          %when there is no noise, they all have equal likelihood
            log_PO_H            = (1/Sim.N)*oney;
        else                                    %otherwise, likelihood is N[C_t, O_t, sig2_o]
            log_PO_H            = -0.5*(S.C(:,t)-R.O(t)).^2/B.sig2_o;
        end
    else                                        %when no observations are present, P[O|H^{(i)}] are all equal
        log_PO_H            = (1/Sim.N)*oney;
    end
    log_n           = oney;                     %extize log sampling spikes
    log_n(sp)       = log(S.p(sp,t));           %compute log P(spike)
    log_n(nosp)     = log(1-S.p(nosp,t));       %compute log P(no spike)

    log_C_Cn        = -0.5*(S.C(:,t)-(B.a*S.C(:,t-1)+B.beta*S.n(:,t))).^2/B.sig2_c;%log P[C_k | C_{t-1}, n_k]

    log_q_n         = oney;                     %initialize log q_n
    log_q_n(sp)     = log(q_n(sp,1));           %compute what was the log prob of sampling a spike
    log_q_n(nosp)   = log(1-q_n(nosp,1));       %or sampling no spike

    log_quotient    = log_PO_H + log_n + log_C_Cn - log_q_n - log_q_C;

    sum_logs        = log_quotient+log(S.w_f(:,t-1));   %update log(weights)
    w               = exp(sum_logs-max(sum_logs));      %exponentiate log(weights)
    S.w_f(:,t)      = w./sum(w);                        %normalize such that they sum to unity

    %% at observations
    if mod(t,Sim.freq)==0
        S = smc_em_bern_stratresamp_v6(Sim,S,t,U_resamp);   %stratified resample
        O = UpdateMoments(Sim,R,B,S,O,t);                   %estimate P[O_s | C_tt] for all t'<tt<s as a gaussian
        ss = t;                                 %store time of last observation
        if mod(t,100)==0 && t<1000              %print # of observations
            fprintf('\b\b\b%d',t)
        elseif mod(t,100)==0 && t<10000
            fprintf('\b\b\b\b%d',t)
        elseif mod(t,100)==0 && t<100000
            fprintf('\b\b\b\b\b%d',t)
        end
    end
end %for time loop

end %function

%% update moments function
function O = UpdateMoments(Sim,R,B,S,O,t)

s               = Sim.freq;                     %find next observation time
O.mu_o(s)       = R.O(t+s);                     %initialize mean of P[O_s | C_s]
O.sig2_o(s)     = B.sig2_o;                     %initialize var of P[O_s | C_s]

if Sim.M>0
    hhat        = zeros(Sim.freq,Sim.M);        %extize hhat
    phat        = zeros(1,Sim.freq+1);          %extize phat

    hs          = S.h(:,t,:);                   %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)= hs(:,1,1:Sim.M);              %this too
    hhat(1,:)   = sum(repmat(S.w_f(:,t),1,Sim.M).*h,1);%initialize hhat
    phat(1)     = sum(S.w_f(:,t).*S.p(:,t),1);  %initialize phat
end

if Sim.M>0
    for tt=1:s
        % update hhat
        for m=1:Sim.M                           %for each spike history term
            hhat(tt+1,m)=B.g(m)*hhat(tt,m)+phat(tt);
        end
        y_t     = B.kx(tt+t)+B.omega'*hhat(tt+1,:)';%input to neuron
        phat(tt+1)  = 1-exp(-exp(y_t)*Sim.dt);      %update phat
    end
else
    phat  = 1-exp(-exp(B.kx(t+1:t+s))*Sim.dt);      %update phat
end

for tt=s:-1:2
    O.mu_n(tt)      = B.beta*phat(:,tt);                    %expected mean of spiking
    O.sig2_n(tt)    = B.beta^2*phat(:,tt).*(1-phat(:,tt));  %expected var of spiking

    O.sig2_c(tt)    = B.sig2_c+O.sig2_n(tt);                %variance of C

    O.mu_o(tt-1)    = B.a^(-1)*(O.mu_o(tt)-O.mu_n(tt));     %mean of P[O_s | C_k]
    O.sig2_o(tt-1)  = B.a^(-2)*(O.sig2_c(tt)+O.sig2_o(tt)); %var of P[O_s | C_k]
end

end %function UpdateMoments