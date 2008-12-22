function [x varargout] = GetDataParams2(varargin)
% this is a hacky function that estimates the "true" parameters of a neuron
% given both a fluorescence trace and the spike times, by optimizing the
% following equations:
% 
% F = alpha*C/(C+k_d)               
% C_t = a C_{t-1} + A n_t + C_0
% 
% note that we assume we know the parameters governing the hill equation
% also note that a = 1-dt/tau
% 
% input: could be either only D, or both D and x0:
%   D: structure with two fields: F and n, both are vectors, could be
%   different lenghts
%   x0: initial parameter guess, x0 = [A a C_0 alpha beta];
% 
% output:
%   x: best guess of parameters
%   Cest: what calcium would be given those parameters
%   Fest: what fluorescence would be given those parameters

D=varargin{1};
if nargin==1
    A   = 25;
    a   = .99;
    C_0 = 0.348;
    alpha = 1;
    beta = 0;
    x0 = [A a C_0 alpha beta];
elseif nargin==2
    x0=varargin{2};
end

T   = min(length(D.F),length(D.n));
F   = D.F(1:T);
n   = D.n(1:T);

k_d = 200;

x = fmincon(@CvsF,x0,[],[],[],[],zeros(size(x0)),inf*ones(size(x0)));


figure(2), clf, plot(F,'k'), hold on, 

Cest = filter(x(1),[1 -x(2)],n)+x(3);
Fest = x(4)*C./(C+k_d) + x(5);
plot(Fest)
stem(n)

function e = CvsF(x)

C = filter(x(1),[1 -x(2)],n)+x(3);
e = norm(F - x(4)*C./(C+k_d) - x(5));

end

end