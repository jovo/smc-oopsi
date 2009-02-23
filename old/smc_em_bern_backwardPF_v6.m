function w_b = smc_em_bern_backwardPF_v6(Sim,S,B)
% this function does the backwards recursion for particle filtering
% assuming spike history terms

w_b   = 1/Sim.N*ones(Sim.N,Sim.T);
ln_Pn = zeros(Sim.N,1);
oney  = ones(Sim.N,1);
zeroy2= zeros(Sim.N);
for t=Sim.T-1:-1:2

    ln_Pn(S.n(:,t)==1)  = log(S.p(S.n(:,t)==1,t));
    ln_Pn(~S.n(:,t))    = log(1-S.p(~S.n(:,t),t));

    C1          = S.C(:,t); 
    C1          = C1(:,oney);
    C0          = B.a*S.C(:,t-1)+B.beta*S.n(:,t);
    C0          = C0(:,oney)';
    ln_PC_Cn    = -0.5*(C0 - C1).^2/B.sig2_c;
        
    ln_Ph_hn = zeroy2;
    for l=1:Sim.M
        h1 = S.h(:,t,l);
        h1 = h1(:,oney);
        h0 = B.g(l)*S.h(:,t-1,l)+S.n(:,t-1);
        h0 = h0(:,oney)';
        ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/B.sig2_h(l);
    end

    sum_lns=ln_Pn(:,oney)+ln_PC_Cn + ln_Ph_hn;
    mx=max(sum_lns,[],1);
    mx=mx(oney,:);
    T0=exp(sum_lns-mx);
    Tn=sum(T0,1);
    T=T0./Tn(oney,:);
    
    PHHn = (T*S.w_f(:,t))';
    PHHn2 = PHHn(oney,:)';
    PHH =  T .* (w_b(:,t)*S.w_f(:,t)')./PHHn2;

    w_b(:,t-1)= sum(PHH,1);
end

end %function