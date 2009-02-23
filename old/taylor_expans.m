clear, clc
C=0:.01:10;
n=1.3;
k_d= 1;
alpha=1;
beta=2;
S=C.^n./(C.^n+k_d);
F=alpha*S+beta;

eps=min(C);
finv=((k_d.*(beta-F))./(F-beta-alpha)).^(1/n)+eps;

% F=1./(1+exp(-2*C+2));
semilogx(C,f,'m')
figure(1), clf, hold on

Finv=((F.*k_d)./(1-F)).^(1/n);

syms k_d beta F alpha n xi
finv=((k_d*(beta-F))/(F-beta-alpha))^(1/n);
dfa=simplify(diff('((k_d*(beta-F))/(F-beta-alpha))^(1/n)',F));
dffa=simplify(diff(dfa,F));
pretty(dfa)
pretty(dffa)


finv=(F*k_d/(1-F))^(1/n);
dfb=simplify(diff('(F*k_d/(1-F))^(1/n)',F));
dffb=simplify(diff(dfb,F));
