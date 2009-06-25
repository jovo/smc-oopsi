function [x C F] = GetDataParams(D)
% this function takes as input a structure D with two fields: F and n, both
% a Tx1 vector. it then fits a model with the following parameters:
% 
% F = alpha*C/(C+k_d)
% C_t = a C_{t-1} + A n_t + C_0

T   = min(length(D.F),length(D.n));
F   = D.F(1:T);
n   = D.n(1:T);

A   = 25;
a   = .99;
C_0 = 0.348;
k_d = 200;
alpha = 1;
beta = 0;

x0 = [A a C_0 alpha beta];
x = fmincon(@CvsF,x0,[],[],[],[],zeros(size(x0)),inf*ones(size(x0)));


C = filter(x(1),[1 -x(2)],n)+x(3);
F = x(4)*C./(C+200)+x(5);

figure(2), clf, plot(D.F(1:T),'k'), hold on, plot(F)

function e = CvsF(x)

C = filter(x(1),[1 -x(2)],n)+x(3);
e = norm(F - x(4)*C./(C+200) - x(5));

end

end