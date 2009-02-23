function w_b = smc_em_bern_backwardPF_v3(Sim,S,B)

% this function does the backwards recursion for particle filtering
% assuming spike history terms

w_b   = 1/Sim.N*ones(Sim.N,Sim.K);
ln_Pn = zeros(Sim.N,1);

for k=Sim.K-1:-1:2

    ln_Pn(S.I(:,k)==1)  = log(S.p(S.I(:,k)==1,k));
    ln_Pn(~S.I(:,k))    = log(1-S.p(~S.I(:,k),k));

    ln_PC_Cn    = -0.5*(meshgrid(S.C(:,k))' - meshgrid(B.a*S.C(:,k-1)+B.beta*S.I(:,k))).^2/B.sig2_c;
    
    ln_Ph_hn = zeros(Sim.N);
    for l=1:Sim.M
        ln_Ph_hn = ln_Ph_hn - 0.5*(meshgrid(S.h(:,k,l))' - meshgrid(B.g(l)*S.h(:,k-1,l)+S.I(:,k-1))).^2/B.sig2_h(l);
    end

    sum_lns=meshgrid(ln_Pn)'+ln_PC_Cn + ln_Ph_hn;    
    mx=max(sum_lns,[],1);
    mx=repmat(mx,Sim.N,1);
    T0=exp(sum_lns-mx);
    T=T0./meshgrid(sum(T0,1));
    
    PHH =  T .* (w_b(:,k)*S.w_f(:,k)')./meshgrid((T*S.w_f(:,k))')';

    w_b(:,k-1)= sum(PHH,1);
    
    if any(isnan(PHH(:)))
        ass=4;
    end
end

end %function