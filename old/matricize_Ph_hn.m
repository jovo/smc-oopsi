ln_Ph_hn = zeroy2;
for m=1:Sim.M
    h1 = S.h(:,t,m);
    h1 = h1(:,oney);
    h0 = B.g(m)*S.h(:,t-1,m)+S.n(:,t-1);
    h0 = h0(:,oney)';
    ln_Ph_hn = ln_Ph_hn - 0.5*(h0 - h1).^2/B.sig2_h(m);
end

mx=max(ln_Ph_hn,[],1);
mx=mx(oney,:);
Ph_hn0=exp(ln_Ph_hn-mx);
Ph_hnn=sum(Ph_hn0,1);
Ph_hn=Ph_hn0./Ph_hnn(oney,:);


ln_Ph_hn2=zeros(Sim.N);
Ph_hn2=zeros(Sim.N);
for i=1:Sim.N
    for j=1:Sim.N
        ln_Ph_hn2(i,j)=...
            -0.5*(S.h(i,t,m)-B.g(m)*S.h(j,t-1,m)-S.n(j,t-1))^2/B.sig2_h(m);
        %         Ph_hn2(i,j)=normpdf(S.h(i,t,m),B.g(m)*S.h(j,t-1,m)-S.n(j,t-1),sqrt(B.sig2_h(m)));
    end
end
norm(ln_Ph_hn2-ln_Ph_hn)
