function R = smc_em_bern_real_exp_v5(Sim,P)
%% this function simulates a neuron according to the following model
% y_t   = k' * x_t + omega' * h_t
% p_t   = p(spike | y_t) = 1-exp(-exp(y_t))
% h_t   = (1-dt/tau_h) h_{t-1} + n_t + epsilon_{ht}, for all h
% C_t   = (1-dt/tau_c) C_{t-1} + beta*I_t + epsilon_{ct}
% F_t   =  alpha*S(C_t) + beta + sigma_ct epsilon_t, epsilon_t ~ N[0,1]
% sig_ct=  gamma*S(C_t) + zeta
% 
% S(C_t) calls Hill_v1 which implements C^n/(C^n+k_d)
% 
% inputs
% Sim   : simulation parameters (eg, dt, Nparticles, etc.)
% P     : parameters of "real" neuron 
% 
% outputs
% R     : states of "real" neuron


%% initialze misc stuff for "real" data
kx        = P.k'*Sim.x;                                 %external input to neuron
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
U_sampl   = rand(1,Sim.T);                              %generate random number to use for sampling
R.C       = P.C_init*ones(1,Sim.T);                     %initialize calcium

if Sim.M>0                                              %if spike history terms, recursively
    R.p         = zeros(1,Sim.T);                       %extize p_t because it must be updated iteratively
    R.n         = zeros(1,Sim.T);                       %extize n_t because it must be updated iteratively
    R.h         = zeros(Sim.M,Sim.T);                   %extize spike history because it must be updated iteratively
    epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(Sim.M,Sim.T); %generate noise on spike history
    for t=2:Sim.T                                       %update states
        R.h(:,t)= (1-Sim.dt./P.tau_h).*R.h(:,t-1)+R.n(t-1) + epsilon_h(:,t);%update h terms
        y_t=kx(t)+P.omega'*R.h(:,t);                    %generate operand for rate function
        R.p(t)=1-exp(-exp(y_t)*Sim.dt);                 %generate rate
        R.n(t)  = U_sampl(t)<R.p(t);                    %sample from bernoulli with prob p_t
    end %time loop
else
    R.p = 1-exp(-exp(kx)*Sim.dt);                       %compute P[n_t]
    R.n = U_sampl<R.p;                                  %sample n_t
end

P.a=Sim.dt/P.tau_c;                                     %for clarity of code
for t=2:Sim.T                                           %recursively update calcium
    R.C(t)  = (1-P.a)*R.C(t-1) + P.A*R.n(t) + P.a*P.C_0 + epsilon_c(t);   
end

F_mu        = P.alpha*Hill_v1(P,R.C)+P.beta;            %compute E[F_t]
F_var       = P.gamma*Hill_v1(P,R.C) + P.zeta;          %compute V[F_t]
R.F         = F_mu+sqrt(F_var).*randn(1,Sim.T);         %add noise to observations
R.F(R.F<0)  = eps;                                      %observations must be non-negative
end