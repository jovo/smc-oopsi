

%% set simulation parameters
Sim.dt      = 0.005;            %time step size (sec)
Sim.T       = 1000;              %# time steps
Sim.x       = ones(1,Sim.T);    %stimulus
Sim.M       = 1;                %# spike history terms

%% set "real" parameters
P.k         = 3;                %bias term
P.tau_c     = 0.5;              %decay rate of calcium
P.beta      = 1;                %jump size of calcium after spike
P.sigma_c   = .1;               %std of noise on calcium
P.sigma_o   = 10*sqrt(P.sigma_c^2*Sim.dt);%std of noise on observations

P.omega     = -1;               %jump size for h after spike
P.tau_h     = 1;                %decay rate for spike history terms
P.sigma_h   = 0.01;             %std of noise on spike history terms

%% initialze misc stuff for "real" data
kx          = P.k'*Sim.x;                               %external input to neuron
epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);    %generate noise on calcium
epsilon_o   = P.sigma_o*randn(1,Sim.T);                 %generate noise on observations
U_sampl     = rand(1,Sim.T);                            %generate random number to use for sampling
R.C         = zeros(1,Sim.T);                           %initialize calcium

if Sim.M>0                                              %if spike history terms, recursively
    R.p         = zeros(1,Sim.T);                       %extize p_t because it must be updated iteratively
    R.n         = zeros(1,Sim.T);                       %extize n_t because it must be updated iteratively
    R.h         = zeros(Sim.M,Sim.T);                   %extize spike history because it must be updated iteratively

    R.p2        = zeros(1,Sim.T);                       %extize p_t because it must be updated iteratively
    R.h2        = zeros(Sim.M,Sim.T);                   %extize spike history because it must be updated iteratively
    epsilon_h   = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(Sim.M,Sim.T); %generate noise on spike history
    for t=2:Sim.T                                       %update states
        R.h(:,t)= (1-Sim.dt./P.tau_h).*R.h(:,t-1)+R.n(t-1) + epsilon_h(:,t);%update h terms
        y_t=kx(t)+P.omega'*R.h(:,t);                    %generate operand for rate function
        R.p(t)=1-exp(-exp(y_t)*Sim.dt);                 %generate rate

        R.h2(:,t)= (1-Sim.dt./P.tau_h).*R.h2(:,t-1)+P.omega'*R.n(t-1) + epsilon_h(:,t);%update h terms
        y_t2=kx(t)+R.h(:,t);                    %generate operand for rate function
        R.p2(t)=1-exp(-exp(y_t2)*Sim.dt);                 %generate rate

        R.n(t)  = U_sampl(t)<R.p(t);                    %sample from bernoulli with prob p_t
    end %time loop
else
    R.p = 1-exp(-exp(kx)*Sim.dt);                      %compute P[n_t]
    R.n = U_sampl<R.p;                                  %sample n_t
end

for t=2:Sim.T                                           %recursively update calcium
    R.C(t)  = (1-Sim.dt/P.tau_c)*R.C(t-1) + P.beta*R.n(t) + epsilon_c(t);
end
R.O = R.C + epsilon_o;                                  %add noise to observations

figure(1), clf, hold on, plot(P.omega*R.h,'k'), plot(R.h2,'r')
norm(P.omega*R.h(1:end-1)-R.h2(1:end-1))