function [n P]=WienerFiltD(F,dt,E)
% this function solves the following optimization problem:
% n = argmin sum_t ((F_t - C_t)^2 + lambda n_t^2
% where
% C_t = a C_{t-1} + n_t, where a=(1-dt/tau)
%
% and then estimates the parameters using:
%
% a = argmin_a ||Wa-Y||^2, where
%   W = -[C_1, C_{T-1}]
%   Y = [{-F_2 n_2}, {-F_T n_T}];
%
% the approach is called ridge regression.  it is essentially the least
% squares solution, regularized for sparsity.
%
% Input----
% F:    fluorescence time series
% dt:   time step size
% E.    structure of parameter estimates
%   tau:time constant
%   lam:prior weight = 1/(rate*A*dt)
%   sig:standard deviation of err

T = length(F);                                  %# frames
o = 1+0*F;                                      %init a unity vector
M = spdiags([-(1 - dt/E.tau)*o o], -1:0,T,T);   %matrix transforming calcium into spikes, ie n=M*C
C = o;                                          %initialize calcium
P = E;
n = M*C;
% E.lam =sqrt(E.lam);

lik     = (F-C)'*(F-C)+E.lam*sum(n.^2);            %initialize Likilihood function
old_lik = inf;
conv    = false;

while conv == false

    D = F-C;
    c = 1/(2*E.sig^2);                          %scale of likelihood term
    g = -2*c*D+2*E.lam*((M*C)'*M)';             %gradient
    H = 2*c*speye(T)+2*E.lam*M'*M;              %Hessian
    C = C-H\g;
    N = M*C;

    old_lik = lik;
    lik     = (F-C)'*(F-C)+E.lam*sum(N.^2);

    if lik-old_lik<=-1e-4
        n = N;
        P = E;
        
        W = -C(1:end-1);
        Y = -F(2:end)+N(2:end);
        a = W'*Y/(W'*W);

        E.tau   = dt/(1-a);
        E.sig   = sqrt((F-C)'*(F-C)/T);
%         E.lam   = T*dt/sum(N);
    else
        conv = true;
    end
end
