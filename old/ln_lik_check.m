
[Sim P]     = InitializeStuff;
P.zeta      = 1e-3;
P.gamma     = 1e-3;
R.C         = 1;
T           = 1000;
epsilon_t   = randn(1,T);
R.C         = 1+randn(1,T);
F_mu        = P.alpha*Hill_v1(P,R.C)+P.beta;            %compute E[F_t]
F_var       = (P.gamma*Hill_v1(P,R.C) + P.zeta);                    %compute V[F_t]
R.F         = F_mu+F_var.*epsilon_t;         %add noise to observations

fake_mu     = linspace(0,1,T);
ln_L        = (-log(F_var)/2-(R.F-fake_mu).^2./F_var);
[foo ind]   = max(ln_L);
[ fake_mu(ind)]

