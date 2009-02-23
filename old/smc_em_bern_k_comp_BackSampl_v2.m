function S = smc_em_bern_k_comp_BackSampl_v2(Sim,R,B)
% the function does the backwards sampling particle filter
% notes: this function assumes spike histories are included.  to turn them
% off, make sure that Sim.M=0 (M is the # of spike history terms).
%
% The backward sampler has standard variance,  as approximated typically.
% Each particle has the SAME backwards sampler, initialized at E[h_t]
% this only does spike history stuff if necessary
%
% the model is F_t = f(C) = alpha C^n/(C^n + k_d) + beta + e_t,
% where e_t ~ N[f(C), gamma*f(C)+zeta]
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
S.p         = zeros(Sim.N,Sim.T);               %extize rate
S.n         = zeros(Sim.N,Sim.T);               %extize spike counts
S.C         = zeros(Sim.N,Sim.T);               %extize calcium
S.w_f       = 1/Sim.N*ones(Sim.N,Sim.T);        %extize forward weights
S.w_b       = 1/Sim.N*ones(Sim.N,Sim.T);        %extize forward weights
S.Neff      = 1/Sim.N*ones(1,Sim.T_o);          %extize N_{eff}

% preprocess stuff for stratified resampling
ints        = linspace(0,1,Sim.N+1);            %generate intervals
diffs       = ints(2)-ints(1);                  %generate interval size
A.U_resamp  = repmat(ints(1:end-1),Sim.T_o,1)+diffs*rand(Sim.T_o,Sim.N); %resampling matrix

% extize misc stuff
A.U_sampl   = rand(Sim.N,Sim.T);                %random samples
A.epsilon_c = sqrt(B.sig2_c)*randn(Sim.N,Sim.T);%generate noise on C

if Sim.M>0                                      %if spike histories
    S.h         = zeros(Sim.N,Sim.T,Sim.M);     %extize spike history terms
    A.epsilon_h   = zeros(Sim.N, Sim.T, Sim.M); %generate noise on h
    for m=1:Sim.M                               %add noise to each h
        A.epsilon_h(:,:,m)   = sqrt(B.sig2_h(m))*randn(Sim.N,Sim.T);
    end
else                                            %if not, comput P[n_t] for all t
    S.p         = repmat(1-exp(-exp(B.kx)*Sim.dt)',1,Sim.N)';
end

if Sim.pf~=0                                    %if doing conditional sampling
    % extize stuff needed for conditional sampling
    A.n_sampl     = rand(Sim.N,Sim.T);          %generate random number to use for sampling n_t
    A.C_sampl     = rand(Sim.N,Sim.T);          %generate random number to use for sampling C_t
    A.oney        = ones(Sim.N,1);              %vector of ones to call for various purposes to speed things up
    A.zeroy       = zeros(Sim.N,1);             %vector of zeros

    % extize stuff needed for REAL backwards sampling
    O.p_o       = zeros(2^(Sim.freq-1),Sim.freq);%extize backwards prob
    O.mu_o      = zeros(2^(Sim.freq-1),Sim.freq);%extize backwards mean
    O.sig2_o    = zeros(1,Sim.freq);            %extize backwards variance

    % extize stuff needed for APPROX backwards sampling
    O.p         = zeros(Sim.freq,Sim.freq);     %extize backwards prob
    O.mu        = zeros(Sim.freq,Sim.freq);     %extize backwards mean
    O.sig2      = zeros(Sim.freq,Sim.freq);     %extize backwards var

    % initialize backwards distributions
    s           = Sim.freq;                     %initialize time of next observation
    O.p_o(1,s)  = 1;                            %initialize P[F_s | C_s]

    [mu1 sig1]  = InitializeLik(B,R.F(s));

    O.mu_o(1,s) = mu1;                       %initialize mean of P[O_s | C_s]
    O.sig2_o(s) = sig1;                     %initialize var of P[O_s | C_s]

    O.p(1,s)    = 1;                            %initialize P[F_s | C_s]
    O.mu(1,s)   = mu1;                       %initialize mean of P[O_s | C_s]
    O.sig2(1,s) = sig1;                     %initialize var of P[O_s | C_s]

    if Sim.freq>1
        for tt=s:-1:2                           %generate spike binary matrix
            A.spikemat(:,tt-1) = repmat([repmat(0,1,2^(s-tt)) repmat(1,1,2^(s-tt))],1,2^(tt-2))';
        end
        nspikes = sum(A.spikemat')';              %count number of spikes at each time step

        for n=0:Sim.freq-1
            A.ninds{n+1}= find(nspikes==n);     %get index for each number of spikes
            A.lenn(n+1) = length(A.ninds{n+1}); %find how many spikes
        end
    end

    O = UpdateMoments_v2(Sim,R,B,S,O,A,s);        %recurse back to get P[O_s | C_s] before the first observation
end

%% loop-de-loop
for t=Sim.freq+1:Sim.T-Sim.freq

    if Sim.pf==0 || R.F(s)>0.98
        S = pf_prior_v2(Sim,R,B,S,A,t);
    else
        S = pf_cond_v2(Sim,R,B,S,O,A,t,s);
    end

    % at observations
    if mod(t,Sim.freq)==0
        S = smc_em_bern_stratresamp_v8(Sim,S,t,A.U_resamp);   %stratified resample
        O = UpdateMoments_v2(Sim,R,B,S,O,A,t);                   %estimate P[O_s | C_tt] for all t'<tt<s as a gaussian
        s = t;                                              %store time of last observation
    end

    if mod(t,100)==0
        if t<1000                          %print # of observations
            fprintf('\b\b\b%d',t)
        elseif t<10000
            fprintf('\b\b\b\b%d',t)
        elseif t<100000
            fprintf('\b\b\b\b\b%d',t)
        end
    end

    if sum(S.C(:,t).^2)==0 || any(isnan(S.C(:,t)))
        ass=4;
    end
end %for time loop

end %conditional sampler

%% initialize likelihood
function [mu1 sig1] = InitializeLik(B,F)
warning off all

% compute mean
finv    = ((B.k_d*(F-B.beta))./(B.alpha-F+B.beta)).^(1/B.n); %initialize search with f^{-1}(o)
mu1     = finv;
% mu1     = fminunc(@fnlogL,finv);        %max P(O|H)
%symbolically compute variance
syms Cest
logL    = -fnlogL(Cest);
dlogL   = diff(logL,'Cest');            %dlog L / dC
ddlogL  = diff(dlogL,'Cest');           %ddlog L / dCC
VC    = -1/ddlogL;                    %neg inverse

% compute variance
Cest    = mu1;                         %eval at max P(O|H)
sig1    = eval(VC);                  %variance approximation

    function logL = fnlogL(C)        %this function compute log L = log P(O|H)
        logL = (((F-fmu_F(C)).^2)./fvar_F(C)+log(fvar_F(C)))/2;
    end

    function mu_F = fmu_F(C)        %this function compute E[F]=f(C)
        mu_F    = B.alpha*C.^B.n./(C.^B.n+B.k_d)+B.beta;
    end

    function var_F = fvar_F(C)      %this function compute V[F]=f(C)
        var_F   = B.gamma*C.^B.n./(C.^B.n+B.k_d)+B.zeta;
    end


end %InitializeLik


%% update moments
function O = UpdateMoments_v2(Sim,R,B,S,O,A,t)

s               = Sim.freq;                     %find next observation time

[mu1 sig1]  = InitializeLik(B,R.F(t+s));

O.mu_o(1,s) = mu1;                       %initialize mean of P[O_s | C_s]
O.sig2_o(s) = sig1;                     %initialize var of P[O_s | C_s]

O.p(1,s)    = 1;                            %initialize P[F_s | C_s]
O.mu(1,s)   = mu1;                       %initialize mean of P[O_s | C_s]
O.sig2(1,s) = sig1;                     %initialize var of P[O_s | C_s]

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
        y_t         = B.kx(tt+t)+B.omega'*hhat(tt+1,:)';%input to neuron
        phat(tt+1)  = 1-exp(-exp(y_t)*Sim.dt);  %update phat
    end
else
    phat  = 1-exp(-exp(B.kx(t+1:t+s)')*Sim.dt); %update phat
end

for tt=s:-1:2
    O.p_o(1:2^(s-tt+1),tt-1)    = repmat(O.p_o(1:2^(s-tt),tt),2,1).*[(1-phat(tt))*ones(1,2^(s-tt)) phat(tt)*ones(1,2^(s-tt))]';
    O.mu_o(1:2^(s-tt+1),tt-1)   = B.a^(-1)*(repmat(O.mu_o(1:2^(s-tt),tt),2,1)-B.beta*A.spikemat(1:2^(s-tt+1),tt-1));     %mean of P[O_s | C_k]
    O.sig2_o(tt-1)              = B.a^(-2)*(B.sig2_c+O.sig2_o(tt)); %var of P[O_s | C_k]

    for n=0:s-tt+1
        nind=A.ninds{n+1};
        O.p(n+1,tt-1)   = sum(O.p_o(nind,tt-1));
        ps              = (O.p_o(nind,tt-1)/O.p(n+1,tt-1))';
        O.mu(n+1,tt-1)  = ps*O.mu_o(nind,tt-1);
        O.sig2(n+1,tt-1)= O.sig2_o(tt-1) + ps*(O.mu_o(nind,tt-1)-repmat(O.mu(n+1,tt-1)',A.lenn(n+1),1)).^2;
    end
end
% sum_n   = A.zeroy;

end %function UpdateMoments


%% this function does particle filtering using the prior sampler
function S = pf_prior_v2(Sim,R,B,S,A,t)

if Sim.M>0                                      %update noise on h
    for m=1:Sim.M
        S.h(:,t,m)=B.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end

    % update rate and sample spikes
    hs              = S.h(:,t,:);               %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_t             = B.kx(t)+B.omega'*h';      %input to neuron
    S.p(:,t)        = 1-exp(-exp(y_t)*Sim.dt);  %update rate for those particles with y_t<0
end
S.n(:,t)            = A.U_sampl(:,t)<S.p(:,t);    %sample

% sample C
S.C(:,t)        = B.a*S.C(:,t-1)+B.A*S.n(:,t)+A.epsilon_c(:,t);

% stratified resample at every observation
if mod(t,Sim.freq)==0
    F_mu        = GenHill(B,S.C(:,t));               %compute E[F_t]
    F_var       = B.gamma*F_mu + B.zeta;          %compute V[F_t]
    ln_w        = -0.5*(R.F(t)-F_mu).^2./F_var + log(F_var)/2;     %compute log of weights
    ln_w        = ln_w-max(ln_w);                       %subtract the max to avoid rounding errors
    w           = exp(ln_w);                            %exponentiate to get actual weights
    S.w_f(:,t)  = w/sum(w);                             %normalize to define a legitimate distribution
    S = smc_em_bern_stratresamp_v8(Sim,S,t,A.U_resamp);
end

end


%% this function does particle filtering using the conditional sampler
function S = pf_cond_v2(Sim,R,B,S,O,A,t,s)

% if spike histories, sample h and update p
if Sim.M>0                                  %update noise on h
    for m=1:Sim.M                           %for each spike history term
        S.h(:,t,m)=B.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end
    hs              = S.h(:,t,:);           %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);      %this too
    S.p(:,t)        = 1-exp(-exp(B.kx(t)+B.omega'*h')*Sim.dt);%update p
end

% compute P[n_k | h_k]
ln_n    = [log(S.p(:,t)) log(1-S.p(:,t))];  %compute [log(spike) log(no spike)]

% compute log G_n(n_k | O_s) for n_k=1 and n_k=0
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

% compute q(n_k | h_k, O_s)
ln_q_n  = ln_n + ln_G;                      %log of sampling dist
mx      = max(ln_q_n,[],2);                 %find max of each column
mx2     = repmat(mx,1,2);                   %matricize
q_n     = exp(ln_q_n-mx2);                  %subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
q_n     = q_n./repmat(sum(q_n,2),1,2);      %normalize to make this a true sampling distribution (ie, sums to 1)

% sample n
S.n(:,t)= A.n_sampl(:,t)<q_n(:,1);          %sample n
% sum_n   = sum_n+S.n(:,t);                   %total number of spikes so far
sp      = S.n(:,t)==1;                      %store index of which samples spiked
nosp    = S.n(:,t)==0;                      %and which did not

% sample C
if mod(t,Sim.freq)==0                       %if not intermittent
    v       = repmat(O.sig2(1,t-s),Sim.N,1);%get var
    m       = repmat(O.mu(1,t-s),Sim.N,1);  %get mean
else
    [fo,sp_i]   = histc(A.C_sampl(sp,t),[0  cumsum(O.p_o(1:k-1,t-s))'/sum(O.p_o(1:k-1,t-s))]);
    [fo,nosp_i] = histc(A.C_sampl(nosp,t),[0  cumsum(O.p_o(1:k,t-s))'/sum(O.p_o(1:k,t-s))]);

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

% update weights
if mod(t,Sim.freq)==0                       %at observations compute P(O|H)
    F_mu        = GenHill_v2(B,S.C(:,t));               %compute E[F_t]
    F_var       = B.gamma*F_mu + B.zeta;          %compute V[F_t]
    log_PO_H    = -0.5*(R.F(t)-F_mu).^2./F_var + log(F_var)/2;     %compute log of weights
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

end %condtional sampler