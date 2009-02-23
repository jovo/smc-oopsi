function O = UpdateMoments(Sim,R,B,S,O,t)
% update moments function

s               = Sim.freq;                     %find next observation time

%% initialize mean and variance
% compute mean
finv            = ((B.k_d*(R.F(t+s)-B.beta))./(B.alpha-R.F(t+s)+B.beta)).^(1/B.n); %initialize search with f^{-1}(o)
ffit            = fminunc(@fnlogL,finv);        %max P(O|H)

% compute variance
syms Cest
logL            = -fnlogL(Cest);
dlogL           = diff(logL,'Cest');            %dlog L / dC
ddlogL          = diff(dlogL,'Cest');           %ddlog L / dCC
VC              = -1/ddlogL;                    %neg inverse
Cest            = ffit;                         %eval at max P(O|H)
varF            = eval(VC);                     %variance approximation

O.mu_o(1,s)     = ffit;                         %initialize mean of P[O_s | C_s]
O.sig2_o(1,s)   = varF;                         %initialize var of P[O_s | C_s]

O.mu(1,s)       = O.mu_o(1,s);                  %initialize mean of P[O_s | C_s]
O.sig2(1,s)     = O.sig2_o(s);                  %initialize var of P[O_s | C_s]

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
    O.mu_o(1:2^(s-tt+1),tt-1)   = B.a^(-1)*(repmat(O.mu_o(1:2^(s-tt),tt),2,1)-B.beta*spikemat(1:2^(s-tt+1),tt-1));     %mean of P[O_s | C_k]
    O.sig2_o(tt-1)              = B.a^(-2)*(B.sig2_c+O.sig2_o(tt)); %var of P[O_s | C_k]

    for n=0:s-tt+1
        nind=A.ninds{n+1};
        O.p(n+1,tt-1)   = sum(O.p_o(nind,tt-1));
        ps              = (O.p_o(nind,tt-1)/O.p(n+1,tt-1))';
        O.mu(n+1,tt-1)  = ps*O.mu_o(nind,tt-1);
        O.sig2(n+1,tt-1)= O.sig2_o(tt-1) + ps*(O.mu_o(nind,tt-1)-repmat(O.mu(n+1,tt-1)',A.lenn(n+1),1)).^2;
    end
end
sum_n   = A.zeroy;

    function logL = fnlogL(C)        %this function compute log L = log P(O|H)
        logL = (((R.F(t+s)-fmu_F(C)).^2)./fvar_F(C)+log(fvar_F(C)))/2;
    end

    function mu_F = fmu_F(C)        %this function compute E[F]=f(C)
        mu_F    = B.alpha*C.^B.n./(C.^B.n+B.k_d)+B.beta; 
    end

    function var_F = fvar_F(C)      %this function compute V[F]=f(C)
        var_F   = B.gamma*C.^B.n./(C.^B.n+B.k_d)+B.zeta;
    end

end %function UpdateMoments
