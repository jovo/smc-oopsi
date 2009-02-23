function w_b = smc_em_bern_backwardPF_v9(Sim,S,B)
% this function does the backwards recursion for particle filtering
% assuming spike history terms

w_b   = 1/Sim.N*ones(Sim.N,Sim.T);
ln_Pn = zeros(Sim.N,1);
oney  = ones(Sim.N,1);
zeroy = zeros(Sim.N);
for t=Sim.T-Sim.freq-1:-1:Sim.freq+1

    ln_Pn(S.n(:,t)==1)  = log(S.p(S.n(:,t)==1,t));
    ln_Pn(~S.n(:,t))    = log(1-S.p(~S.n(:,t),t));

    C1          = S.C(:,t);
    C1          = C1(:,oney);
    C0          = (1-B.a)*S.C(:,t-1)+B.Aeff*S.n(:,t)+B.a*B.C_0;
    C0          = C0(:,oney)';
    ln_PC_Cn    = -0.5*(C0 - C1).^2/B.sig2_c;

    ln_Ph_hn = zeroy;
    for m=1:Sim.M
        h1 = S.h(:,t,m);
        h1 = h1(:,oney);
        h0 = B.g(m)*S.h(:,t-1,m)+S.n(:,t-1);
        h0 = h0(:,oney)';
        ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/B.sig2_h(m);
    end

    sum_lns = ln_Pn(:,oney)+ln_PC_Cn + ln_Ph_hn;
    mx      = max(sum_lns,[],1);
    mx      = mx(oney,:);
    T0      = exp(sum_lns-mx);
    Tn      = sum(T0,1);
    T       = T0./Tn(oney,:);

    PHHn    = (T*S.w_f(:,t-1))';
    PHHn(PHHn==0) = eps;
    PHHn2   = PHHn(oney,:)';
    PHH     =  T .* (w_b(:,t)*S.w_f(:,t-1)')./PHHn2;
    PHH     =  PHH/sum(PHH(:));

    w_b(:,t-1)= sum(PHH,1);

    if mod(t,100)==0 && t>=9900
        fprintf('\b\b\b\b\b%d',t)
    elseif mod(t,100)==0 && t>=900
        fprintf('\b\b\b\b%d',t)
    elseif mod(t,100)==0
        fprintf('\b\b\b%d',t)
    end
end

end %function