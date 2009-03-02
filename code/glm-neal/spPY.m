function Z=spPY(P,Y,X,t)
%HMM observation function P(Y|X,t).
%Z=PY(p,[],X,t) computes a sample Y from P(Y|X,t). If X is a NxM matrix of 
%M dimensional states, Z is a NxK vector of the samples.
%Z=PY(p,Y,X,t) computes P(Y|X,t). If X is a NxM matrix of N M-dimensional 
%states, Z is a Nx1 vector of the probabilities. "p" is a structure of 
%parameters. YM(01/2009)

%REQUIRED PARAMETERS
% P.alpha=1;    %fluorescence scale
% P.beta=0;     %fluorescence offset
% P.gamma=1e-4; %fluorescence conversion factor
% P.zeta=4e-4;  %fluorescence background noise
% P.n=1;        %fluorescence saturation exponent
% P.k_d=100;    %fluorescence saturation offset

if(isempty(Y))  
  C=X(:,2);                       %[Ca]
  S=Hill_v1(P,C);
  s2=P.gamma*S+P.zeta;            %variance
  S=P.alpha*S+P.beta;             %mean  
  Z=S+sqrt(s2).*randn(size(S));   %observation
else
  C=X(:,2);                       %[Ca]
  PG=X(:,3);                      %PG probs
  S=Hill_v1(P,C);
  s2=P.gamma*S+P.zeta;            %variance, expected
  S=P.alpha*S+P.beta;             %mean, expected  
  dY=(Y-S).^2;                    %difference
  adj=0;
  Z=exp((adj-dY)./(2*s2))./sqrt(2*pi*s2)./PG;
end


function F = Hill_v1(P,C)
% generalized Hill model
C(C<0) = 0; F = C.^P.n./(C.^P.n+P.k_d^P.n);