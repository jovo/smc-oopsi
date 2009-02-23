% clear C0 C1
ln_Pn(S.n(:,t)==1)  = log(S.p(S.n(:,t)==1,t));
ln_Pn(~S.n(:,t))    = log(1-S.p(~S.n(:,t),t));

C1          = S.C(:,t);
C1          = C1(:,oney);
C0          = B.a*S.C(:,t-1)+B.beta*S.n(:,t);
C0          = C0(:,oney)';
ln_PC_Cn    = -0.5*(C0 - C1).^2/B.sig2_c;

ln_Ph_hn = zeroy2;
for m=1:Sim.M
    h1 = S.h(:,t,m);
    h1 = h1(:,oney);
    h0 = B.g(m)*S.h(:,t-1,m)+S.n(:,t-1);
    h0 = h0(:,oney)';
    ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/B.sig2_h(m);
end

% mx=max(ln_Ph_hn,[],1);
% mx=mx(oney,:);
% Ph_hn0=exp(ln_Ph_hn-mx);
% Ph_hnn=sum(Ph_hn0,1);
% Ph_hn=Ph_hn0./Ph_hnn(oney,:);
% 
sum_lns=ln_Pn(:,oney)+ln_PC_Cn + ln_Ph_hn;
mx=max(sum_lns,[],1);
mx=mx(oney,:);
T0=exp(sum_lns-mx);
Tn=sum(T0,1);
T=T0./Tn(oney,:);

PHHn = (T*S.w_f(:,t-1))';
PHHn2 = PHHn(oney,:)';
PHH =  T .* (w_b(:,t)*S.w_f(:,t-1)')./PHHn2;

w_b(:,t-1)= sum(PHH,1);

any(isnan(w_b(:,t-1)))