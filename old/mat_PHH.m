PHHn = (T*S.w_f(:,t-1))';
PHHn2 = PHHn(oney,:)';
PHH =  T .* (w_b(:,t)*S.w_f(:,t-1)')./PHHn2;


PHHnb=(T0*S.w_f(:,t-1))';
find(PHHnb==0);

PHHnc=(exp(sum_lns)*S.w_f(:,t-1))';
find(PHHnc==0);

PHHn2b = PHHnb(oney,:)';
PHHb =  T .* (w_b(:,t)*S.w_f(:,t-1)')./PHHn2b;
[foo ind]=find(~isfinite(PHHb))


ass=sum_lns+repmat(log(S.w_f(:,t-1)'),100,1);
norm(exp(sum(ass,2))-exp(sum_lns*S.w_f(:,t-1)))

size(-exp(sum_lns*S.w_f(:,t-1)))

size(exp(sum(ass,1)))



T2=zeros(Sim.N);
ln_Pns2=zeros(Sim.N);
ln_PC_Cn2=zeros(Sim.N);
ln_Ph_hn2=zeros(Sim.N);
for i=1:Sim.N
    for j=1:Sim.N
        T2(i,j)=ln_Pn(i)...
            -0.5*(S.C(i,t)-(B.a*S.C(j,t-1)+B.beta*S.n(i,t)))^2/B.sig2_c...
            -0.5*(S.h(i,t,1)-(B.g(1)*S.h(j,t-1,1)+S.n(j,t-1)))^2/B.sig2_h;

        ln_Pns2(i,j)=ln_Pn(i);
        ln_PC_Cn2(i,j) = -0.5*(S.C(i,t)-(B.a*S.C(j,t-1)+B.beta*S.n(i,t)))^2/B.sig2_c;
        ln_Ph_hn2=-0.5*(S.h(i,t,m)-(B.g(1)*S.h(j,t-1,m)+S.n(j,t-1)))^2/B.sig2_h(m);
    end
end

norm(T-T2)
norm(ln_Pns2-ln_Pn(:,oney))
norm(ln_PC_Cn-ln_PC_Cn2)
norm(ln_Ph_hn-ln_Ph_hn2)

ln_Ph_hn3=- 0.5*(meshgrid(S.h(:,t,1))' - meshgrid(B.g(1)*S.h(:,t-1,1)+S.n(:,t-1))).^2/B.sig2_h(1);
norm(ln_Ph_hn2-ln_Ph_hn3)

h1 = S.h(:,t,m);
h1 = h1(:,oney);
h0a = B.g(m)*S.h(:,t-1,m);
h0a = h0a(:,oney)';

h0b = S.n(:,t-1);
h0b = h0b(:,oney); size(h0b), sum(h0b(:))
ln_Ph_hn4 =  - 0.5*(h0a+h0b-h1).^2/B.sig2_h(m);
norm(ln_Ph_hn2-ln_Ph_hn4)


