function w_b = smc_em_bern_backwardPF_v5(Sim,S,B)
% this function does the backwards recursion for particle filtering
% assuming spike history terms

w_b   = 1/Sim.N*ones(Sim.N,Sim.K);
ln_Pn = zeros(Sim.N,1);
oney  = ones(Sim.N,1);

for k=Sim.K-1:-1:2

    ln_Pn(S.I(:,k)==1)  = log(S.p(S.I(:,k)==1,k));
    ln_Pn(~S.I(:,k))    = log(1-S.p(~S.I(:,k),k));

    C1          = S.C(:,k); 
    C1          = C1(:,oney);
    C0          = B.a*S.C(:,k-1)+B.beta*S.I(:,k);
    C0          = C0(:,oney)';
    ln_PC_Cn    = -0.5*(C0 - C1).^2/B.sig2_c;
        
    ln_Ph_hn = zeros(Sim.N);
    ln_Ph_hn2 = zeros(Sim.N);
    for l=1:Sim.M
        h1 = S.h(:,k,l);
        h1 = h1(:,oney);
        h0 = B.g(l)*S.h(:,k-1,l)+S.I(:,k-1);
        h0 = h0(:,oney)';
        ln_Ph_hn2 = ln_Ph_hn2 - 0.5*(h0 - h1).^2/B.sig2_h(l);
    end

    sum_lns=ln_Pn(:,oney)+ln_PC_Cn + ln_Ph_hn;
    mx=max(sum_lns,[],1);
    mx=mx(oney,:);
    T0=exp(sum_lns-mx);
    Tn=sum(T0,1);
    T=T0./Tn(oney,:);
    
    PHHn = (T*S.w_f(:,k))';
    PHHn2 = PHHn(oney,:)';
    PHH =  T .* (w_b(:,k)*S.w_f(:,k)')./PHHn2;

    w_b(:,k-1)= sum(PHH,1);
end

end %function