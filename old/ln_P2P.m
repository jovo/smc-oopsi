function P = ln_P2P(ln_P)

mx=max(ln_P,[],2);
mx=repmat(mx,1,length(ln_P));
P=exp(ln_P-mx);
P=P./meshgrid(sum(P,1));
