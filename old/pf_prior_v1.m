function S = pf_prior_v1(Sim,R,B,S,A,t)
% this function does particle filtering using the prior sampler

if Sim.M>0                                      %update noise on h
    for m=1:Sim.M
        S.h(:,t,m)=B.g(m)*S.h(:,t-1,m)+S.n(:,t-1)+A.epsilon_h(:,t,m);
    end

    % update rate and sample spikes
    hs              = S.h(:,t,:);               %this is required for matlab to handle a m-by-n-by-p matrix
    h(:,1:Sim.M)    = hs(:,1,1:Sim.M);          %this too
    y_t             = B.kx(t)+B.omega'*h';      %input to neuron
    S.p(:,t)        = 1-exp(-exp(y_t)*Sim.dt);  %update rate for those particles with y_t<0
end
S.n(:,t)            = A.U_sampl(:,t)<S.p(:,t);    %sample

% sample C
S.C(:,t)        = B.a*S.C(:,t-1)+B.A*S.n(:,t)+A.epsilon_c(:,t);

% stratified resample at every observation
if mod(t,Sim.freq)==0
    F_mu        = GenHill(B,S.C(:,t));               %compute E[F_t]
    F_mu(imag(F_mu)~=0) = 0;
    F_var       = B.gamma*F_mu + B.zeta;          %compute V[F_t]
    ln_w        = -0.5*(R.F(t)-F_mu).^2./F_var + log(F_var)/2;     %compute log of weights
    ln_w        = ln_w-max(ln_w);                       %subtract the max to avoid rounding errors
    w           = exp(ln_w);                            %exponentiate to get actual weights
    S.w_f(:,t)  = w/sum(w);                             %normalize to define a legitimate distribution
    S = smc_em_bern_stratresamp_v8(Sim,S,t,A.U_resamp);
end

end