function S = pf_cond_v1(Sim,R,B,S,O,A,t,s)
% this function does particle filtering using the conditional sampler

%% if spike histories, sample h and update p
if Sim.M>0                                  %update noise on h
    for m=1:Sim.M                           %for each spike history term
        S.h(:,t,m)=B.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end
    hs              = S.h(:,t,:);           %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);      %this too
    S.p(:,t)        = 1-exp(-exp(B.kx(t)+B.omega'*h')*Sim.dt);%update p
end

%% compute P[n_k | h_k]
ln_n    = [log(S.p(:,t)) log(1-S.p(:,t))];  %compute [log(spike) log(no spike)]

%% compute log G_n(n_k | O_s) for n_k=1 and n_k=0
k   = Sim.freq-(t-s)+1;

m0  = B.a*S.C(:,t-1);                       %mean of P[C_k | C_{t-1}, n_k=0]
m1  = B.a*S.C(:,t-1)+B.A;                   %mean of P[C_k | C_{t-1}, n_k=1]

m2  = O.mu(1:k,t-s);                        %mean of P[O_s | C_k] for n_k=1 and n_k=0
v2  = O.sig2(1:k,t-s);                      %var of P[O_s | C_k] for n_k=1 and n_k=0
v   = repmat(B.sig2_c+v2',Sim.N,1);         %var of G_n(n_k | O_s) for n_k=1 and n_k=0

ln_G0= -0.5*log(2*pi.*v)-.5*(repmat(m0,1,k)-repmat(m2',Sim.N,1)).^2./v;   %log G_n(n_k | O_s) for n_k=1 and n_k=0
ln_G1= -0.5*log(2*pi.*v)-.5*(repmat(m1,1,k)-repmat(m2',Sim.N,1)).^2./v;   %log G_n(n_k | O_s) for n_k=1 and n_k=0

mx  = max(max(ln_G0,[],2),max(ln_G1,[],2))';%get max of these
mx  = repmat(mx,k,1)';

G1  = exp(ln_G1-mx);                        %norm dist'n for n=1;
M1  = G1*O.p(1:k,t-s);                      %times prob of n=1

G0  = exp(ln_G0-mx);                        %norm dist'n for n=0;
M0  = G0*O.p(1:k,t-s);                      %times prob n=0

ln_G    = [log(M1) log(M0)];                %ok, now we actually have the gaussians

%% compute q(n_k | h_k, O_s)
ln_q_n  = ln_n + ln_G;                      %log of sampling dist
mx      = max(ln_q_n,[],2);                 %find max of each column
mx2     = repmat(mx,1,2);                   %matricize
q_n     = exp(ln_q_n-mx2);                  %subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
q_n     = q_n./repmat(sum(q_n,2),1,2);      %normalize to make this a true sampling distribution (ie, sums to 1)

%% sample n
S.n(:,t)= A.n_sampl(:,t)<q_n(:,1);          %sample n
sum_n   = sum_n+S.n(:,t);                   %total number of spikes so far
sp      = S.n(:,t)==1;                      %store index of which samples spiked
nosp    = S.n(:,t)==0;                      %and which did not

%% sample C
if mod(t,Sim.freq)==0                       %if not intermittent
    v       = repmat(O.sig2(1,t-s),Sim.N,1);%get var
    m       = repmat(O.mu(1,t-s),Sim.N,1);  %get mean
else
    [fo,sp_i]   = histc(C_sampl(sp,t),[0  cumsum(O.p_o(1:k-1,t-s))'/sum(O.p_o(1:k-1,t-s))]);
    [fo,nosp_i] = histc(C_sampl(nosp,t),[0  cumsum(O.p_o(1:k,t-s))'/sum(O.p_o(1:k,t-s))]);

    v       = O.sig2(1:k,t-s);              %get var
    v(sp)   = v(sp_i);                      %if particle spiked, then use variance of spiking
    v(nosp) = v(nosp_i);                    %o.w., use non-spiking variance

    m       = O.mu(1:k,t-s);                %get mean
    m(sp)   = m(sp_i);                      %if particle spiked, use spiking mean
    m(nosp) = m(nosp_i);                    %o.w., use non-spiking mean
end
v_c         = (1./v+1/B.sig2_c).^(-1);      %variance of dist'n for sampling C
m_c         = v_c.*(m./v+(B.a*S.C(:,t-1)+B.A*S.n(:,t))/B.sig2_c);%mean of dist'n for sampling C
S.C(:,t)    = normrnd(m_c,sqrt(v_c));       %sample C

%% update weights
if mod(t,Sim.freq)==0                       %at observations compute P(O|H)
    log_PO_H            = -0.5*(S.C(:,t)-R.O(t)).^2/B.sig2_o;
else                                        %when no observations are present, P[O|H^{(i)}] are all equal
    log_PO_H            = (1/Sim.N)*A.oney;
end
log_n           = A.oney;                   %extize log sampling spikes
log_n(sp)       = log(S.p(sp,t));           %compute log P(spike)
log_n(nosp)     = log(1-S.p(nosp,t));       %compute log P(no spike)
log_C_Cn        = -0.5*(S.C(:,t)-(B.a*S.C(:,t-1)+B.A*S.n(:,t))).^2/B.sig2_c;%log P[C_k | C_{t-1}, n_k]

log_q_n         = A.oney;                   %initialize log q_n
log_q_n(sp)     = log(q_n(sp,1));           %compute what was the log prob of sampling a spike
log_q_n(nosp)   = log(1-q_n(nosp,1));       %or sampling no spike
log_q_C         = -0.5*(S.C(:,t)-m_c).^2./v_c;  %log prov of sampling the C_k that was sampled

log_quotient    = log_PO_H + log_n + log_C_Cn - log_q_n - log_q_C;

sum_logs        = log_quotient+log(S.w_f(:,t-1));   %update log(weights)
w               = exp(sum_logs-max(sum_logs));      %exponentiate log(weights)
S.w_f(:,t)      = w./sum(w);                        %normalize such that they sum to unity

%% at observations
if mod(t,Sim.freq)==0
    S = smc_em_bern_stratresamp_v8(Sim,S,t,U_resamp);   %stratified resample
    O = UpdateMoments_v2(Sim,R,B,S,O,t);                   %estimate P[O_s | C_tt] for all t'<tt<s as a gaussian
    s = t;                                              %store time of last observation
end

end