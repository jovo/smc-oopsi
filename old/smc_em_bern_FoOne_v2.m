function S = smc_em_bern_FoOne_v2(Sim,R,B)

%%%%%%%% initialize stuff
S.p     = zeros(Sim.N,Sim.K);           %initialize rate
S.I     = zeros(Sim.N,Sim.K);           %initialize spike counts
S.C     = R.C(1)*ones(Sim.N,Sim.K);     %initialize calcium
S.Neff  = 1/Sim.N*ones(1,Sim.K_o);      %initialize N_{eff}

U_sampl     = rand(Sim.N,Sim.K);        %random samples
epsilon_c   = sqrt(B.sig2_c)*randn(Sim.N,Sim.K); %generate noise on c
xk          = B.k'*Sim.x;               %compute input to neuron
sig2_z      = (1/B.sig2_o + 1/B.sig2_c)^(-1);       %var of sample distribution

%for testing purposes only....take out after
ints        = linspace(0,1,Sim.N+1);
diffs       = ints(2)-ints(1);
U_resamp    = repmat(ints(1:end-1),Sim.K_o,1)+diffs*rand(Sim.K_o,Sim.N);
Nresamp     = 0;

    
if Sim.M>0 %if spike history terms are present
    %generate noise on h
    S.h         = R.h(1)*ones(Sim.N,Sim.K,Sim.M);  %initialize spike history terms
    epsilon_h   = zeros(Sim.N, Sim.K, Sim.M);
    for l=1:Sim.M
        epsilon_h(:,:,l)   = sqrt(B.sig2_h(l))*randn(Sim.N,Sim.K);   %generate noise on h
    end

    %preprocess stuff for stratified resampling
    ints        = linspace(0,1,Sim.N+1);
    diffs       = ints(2)-ints(1);
    U_resamp    = repmat(ints(1:end-1),Sim.K_o,1)+diffs*rand(Sim.K_o,Sim.N);
    Nresamp     = 0;
else %if no spike histories, generate p and sample I
    S.p(1,:) = 1-exp(-exp(xk));                                  %update rate for those particles with y_t<0
    S.p      = repmat(S.p(1,:),Sim.N,1);%make rate the same for each particle
    S.I      = U_sampl<S.p;             %sample spikes
end

for k=2:Sim.K

    %%%%%%%% if h's, update h's, update p_k, and sample I_k
    if Sim.M>0

        %update h
        for l=1:Sim.M
            S.h(:,k,l)=B.g(l)*S.h(:,k-1,l)+S.I(:,k-1)+epsilon_h(:,k,l);
        end

        %generate input
        hs              = S.h(:,k,:);               %this is required for matlab to handle a m-by-n-by-p matrix
        h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
        y_t             = xk(k)+B.omega'*h';        %input to neuron
        S.p(:,k)        = 1-exp(-exp(y_t));         %update rate for those particles with y_t<0
    end

    %generate sampling distribution
    if mod(k*Sim.dt,Sim.dtSample)==0                    %at observation times
        ln_I    = [log(S.p(:,k)) log(1-S.p(:,k))];          %compute [log(spike) log(no spike)]
        mu_ck   = [B.a*S.C(:,k-1)+B.beta B.a*S.C(:,k-1)];   %mean when spiked and not spiked

        mu_zk   = sig2_z*(mu_ck/B.sig2_c + R.O(k)/B.sig2_o);  %mean of sample distribution
        ln_Z    = 0.5* (mu_zk.^2/sig2_z - mu_ck.^2/B.sig2_c - R.O(k).^2/B.sig2_o); %log of normalizing constant (not including constants)
        ln_Opt  = ln_I + ln_Z;                              %log of sampling dist
        mx      = max(ln_Opt,[],2);                         %find max of each column
        mx2     = repmat(mx,1,2);                           %matricize
        OptP    = exp(ln_Opt-mx2);                          %subtract max to ensure that for each column, there is at least one positive probability, and exponentiate
        OptP    = OptP./repmat(sum(OptP,2),1,2);            %normalize to make this a true sampling distribution (ie, sums to 1)
        S.I(:,k)= U_sampl(:,k)<OptP(:,1);                   %sample p(spike | y_t)
    else                                                %at non-observation times
        OptP    = 1-S.p(:,k);                           %again, prob of NOT spiking
        S.I(:,k)= U_sampl(:,k)>OptP;                    %sample (note the > because we computed p(no spike | y_t)
    end

    %%%%%%%% sample C
    if mod(k*Sim.dt,Sim.dtSample)==0  %at observation times
        mu_ck   = B.a*S.C(:,k-1)+B.beta*S.I(:,k);   %mean when spiked and not spiked
        mu_zk   = sig2_z*(mu_ck/B.sig2_c + R.O(k)/B.sig2_o);  %mean of sample distribution
        S.C(:,k)=normrnd(mu_zk,sqrt(sig2_z));
    else                                %at non-observation times
        S.C(:,k)=B.a*S.C(:,k-1)+B.beta*S.I(:,k)+epsilon_c(:,k);
    end

    %%%%%%%% stratified resample
        if mod(k*Sim.dt,Sim.dtSample)==0 %&& Sim.M>0 
            S = smc_em_bern_stratresamp_v1(Sim,S,R,k,B.sig2_o,U_resamp);
        end %resample
end %for k=2:Sim.K

end %function S = smc_em_bern_FoOne_v1(Sim,R,B)
