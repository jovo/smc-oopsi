clear; clc

Xim.dt = 1/30;
Xim.T   = 500;
Xim.MaxIter = 10;

PQ.gam  = (1-Xim.dt/1);
H.n    = 1; 
H.k_d  = 200;
PQ.A    = 50;
PQ.zeta = 0.02;

n = rand(1,Xim.T)>.98;
C = filter(1,[1 -PQ.gam],n*PQ.A);               % calcium concentration
F = Hill_v1(H,C)+0.1+PQ.zeta*randn(1,Xim.T);
% F = C+PQ.zeta*randn(1,Xim.T);
figure(1), clf, 
subplot(311), plot(F), axis('tight') 
subplot(312), plot(C), axis('tight') 
subplot(313), bar(n), axis('tight') 

%%
[S.n P2]= fast_oopsi(F,Xim,PQ);
S.nnorm = S.n/max(S.n);
figure(1), subplot(313), hold all, stem(S.nnorm)
S.C     = filter(1,[1 -PQ.gam],S.nnorm'*PQ.A);               % calcium concentration
figure(1), subplot(312), hold all, plot(S.C,'--')
C1      = [Hill_v1(H,S.C); ones(1,Xim.T)];
ab      = C1'\F';
P.alpha = ab(1);
P.beta  = ab(2);
P.zeta  = sqrt(sum((F-ab'*C1).^2)/Xim.T);

[ab; P.zeta/PQ.zeta]
% resid   = [F-ab'*C1; ones(1,Xim.T)];
% gz      = resid'\F';
% P.gamma = gz(1);
% P.zeta  = gz(2);

