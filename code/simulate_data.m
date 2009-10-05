clear; clc

V.dt = 1/30;
V.T   = 500;

P.tau_c = 1;
P.gam   = 1-V.dt/P.tau_c;
P.n     = 1; 
P.k_d   = 200;
P.A     = 50;
P.sig   = 0.02;
P.a     = 1;
P.b     = 0.1;

n = rand(1,V.T)>.98;
C = filter(1,[1 -P.gam],n*P.A);               % calcium concentration
F = P.a*Hill_v1(P,C)+P.b+P.sig*randn(1,V.T);
figure(1), clf, 
subplot(311), plot(F), axis('tight') 
subplot(312), plot(C), axis('tight') 
subplot(313), bar(n), axis('tight') 
