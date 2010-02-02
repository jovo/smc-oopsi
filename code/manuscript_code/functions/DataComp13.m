function I = DataComp13(F,P,Sim)
% this function infers the spike train from the fluorescent time-series
% using one of several algorithms we've developed
% 
% Input---
% F:    vector of fluorescence traces, linearly normalized between 0+eps and 1
% P:    initial parameter estimates
% Sim:  structure containing other info (some shown below)
%   dt:     time step size
%   T:      # time steps
%   Alg:    which alg to use
%   x:      stimulus (optional)
%   freq:   intermittency
%   MaxIter:# EM iterations before stopping
% 
% Output is structure I with the following fields
% name:     name of algorithm used to infer spike train
% n:        inferred spike train (if only mean is computed)
% P:        MLE of parameters
% M:        moments and percentiles when smc-em is used (when using smc-em)

switch Sim.Alg
    case 1 %quadratic optimization
        o           = [1; ones(Sim.T-1,1)];
        M           = spdiags([-(1 - Sim.dt/P.tau)*o o], -1:0,Sim.T,Sim.T);  %matrix transforming calcium into spikes, ie n=M*C
        I.n         = M*F;
        I.P         = P;
        I.name      = [{'Least'};{'Squares'}];
    case 2 %wiener filter
        [I.n I.P]   = WienerFiltD(F,Sim.dt,P);
        I.name      = [{'Wiener Filter'}];
    case 3 %YF06
        [n,f,f]     = CaDeconvStandAlone(F,P.tau,0,0,0,0,5000,Sim.dt,0.38);
        I.name      = [{'Yaksi and'};{'Friedrich'}];
        I.n         = n;
    case 4 %barrier method
        I.name      = [{'Non-negative'};{'Sparse Filter'}];
        [I.n I.P]   = BarrierFiltD(F,P,Sim);
    case 5 %ppr (seems to have a bug)
        I.name      = [{'Projection'};{'Pursuit'};{'Regression'}];
        [I.n I.P]   = PPRFiltD(F,P,Sim);
    case 6 %regularized ppr (seems to have a bug)
        I.name      = [{'Regularized'};{'PPR'}];
        [I.n I.P]   = PPRRegFiltD(F,P,Sim);
    case 7 %particle smoother without stimulus
        I.name      = [{'OOPSI Filter'}];
%         Sim.x       = ones(size(Sim.x));
        [I.M I.P]   = ParticleFilterD(F,P,Sim);
%     case 8 %intermittent
%         I.name      = [{'Intermittent'};{'Particle'};{'Filter'}];
%         freq        = 2;
%         F           = NaN*zeros(1,T*freq);
%         F(freq:freq:end)= F;
%         n_est       = SubSampleSpikes(D,freq);
%         F           = F(1:length(n_est));
%         [I.M I.P]   = ParticleFilterD(F,P,Sim);
%     case 9 %particle filter with stimulus
%         I.name      = [{'Stimulus Conditioned'};{'Particle Smoother'}];
%         [I.M I.P]   = ParticleFilterD(F,P,Sim);
end