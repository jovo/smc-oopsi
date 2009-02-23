function [R P] = smc_em_bern_real_expquad_v1(Sim)
% this function simulates a neuron with spike history terms where
% f(s)=exp(s)
% 
%%%%%%%%%% inputs
% Sim   : simulation parameters (eg, dt, Nparticles, etc.)
% 
%%%%%%%%%% outputs
% R     : states of "real" neuron
% P     : parameters of "real" neuron 

%%%%%%%% set "real" parameters
P.k         = -1*ones(Sim.KernelSize,1);%3*sin(linspace(pi/2,8*pi/2,Sim.KernelSize))'; %bias and stimulus kernel
P.k(1)      = -3;                                           %bias term

P.tau_c     = 0.5;                                          %decay rate of calcium
P.beta      = 2;                                            %jump size of calcium after spike
P.sigma_c   = 1;                                            %std of noise on calcium
P.sigma_o   = 2*sqrt(P.sigma_c^2*Sim.dt);                   %std of noise on observations

if Sim.M>0                                                  %if spike history terms are included
    P.omega     = -0.01*ones(Sim.M,1)-0.001*rand(Sim.M,1);  %jump size for h after spike
    P.tau_h     = 0.2*ones(Sim.M,1)+0.001*rand(Sim.M,1);    %decay rate for spike history terms
    P.sigma_h   = 0.01*ones(Sim.M,1);                       %std of noise on spike history terms
end

%%%%%%%% initialize "real" data
R.p         = zeros(1,Sim.K);       %spike rate
R.I         = zeros(1,Sim.K);       %spike times
R.C         = 1*ones(1,Sim.K);      %initialize calcium
R.O         = NaN*zeros(1,Sim.K);   %initialize observations
if Sim.M>0
    R.h = 20*ones(Sim.M,Sim.K); %spike history
end

%%%%%%%% simulate "real" data
xk        = P.k'*Sim.x;                                 %external input to neuron
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.K);      %generate noise on calcium
epsilon_o = P.sigma_o*randn(1,Sim.K);                   %generate noise on observations
U_sampl   = rand(1,Sim.K);                              %generate random number to use for sampling

if Sim.M>0
    epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.K).*randn(Sim.M,Sim.K); %generate noise on spike history
    for k=2:Sim.K                                       %update states
        R.h(:,k)= (1-Sim.dt./P.tau_h).*R.h(:,k-1)+R.I(k-1) + epsilon_h(:,k);%update h terms
        y_t=xk(k)+P.omega'*R.h(:,k);                    %generate operand for rate function
        R.p(k)=1-exp(-exp(y_t));                        %generate rate
        R.I(k)  = U_sampl(k)<R.p(k);              %sample from poisson with rate proportional to lambda(k)
    end %time loop
else
    R.p = 1-exp(-exp(xk));
    R.I = U_sampl<R.p;
end

for k=2:Sim.K     %update calcium
    R.C(k)  = (1-Sim.dt/P.tau_c)*R.C(k-1) + P.beta*R.I(k) + epsilon_c(k);   
end
obs     = [Sim.dtSample/Sim.dt:Sim.dtSample/Sim.dt:Sim.K];  %vector of observation times
R.O(obs)= R.C(obs) + epsilon_o(obs);                        %add noise to observations

end