function Z=spPXX(P,X1,X0,t)
%HMM transition function P(X1|X0,t). 
%Z=PXX(p,[]) computes a 1xM-dimensional sample from P(X).
%Z=PXX(p,X1) computes the probability P(X1). If X1 is NxM matrix of N 
%M-dimensional states, Z is a Nx1 vector.
%Z=PXX(p,[],X0,t) computes a 1xM-dimensional sample from P(X1|X0,t). 
%Z=PXX(p,X1,X0,t) computes P(X1|X0,t). If X1 and X0 are NxM and 
%KxM matrices of N and K M-dimensional states, Z is a NxK matrix 
%P(i,j)=P(X1(i)|X0(j),t). "p" is a structure of parameters. YM(01/2009)

%REQUIRED PARAMETERS
% P.dt=0.02;        %time bin
% P.k=1;            %spike generation probabilities, may be 1xT array
% P.A=50;           %Ca jump 
% P.sigma_c=15;     %Ca noise
% P.tau_c=0.3;      %Ca reporter time constant
% P.C_0=50;         %Ca background
%[optional]
% P.sig2=P.sigma_c^2*P.dt    %adjusted variance
% P.rate=exp(-P.k*P.dt)     %spiking rate, may be 1xT array
if(nargin<4) t=1; end

if(isfield(P,'sig2')) s2=P.sig2; else s2=P.sigma_c^2*P.dt; end  %precompute
if(isfield(P,'rate'))                                           %precompute
  if(length(P.rate)>=t) p=P.rate(t); else p=P.rate(end); end
else
  if(length(P.k)>=t) p=exp(-P.k(t)*P.dt); else p=exp(-P.k(end)*P.dt); end
end

if(nargin<3 && isempty(X1))
  n=rand>p;                   %spike
  C=P.C_0+sqrt(s2)*randn;     %base Ca
  Z=[n C];
elseif(nargin<3)
  n=X1(:,1); C=X1(:,2);       %unwrap state
  dC=(C-P.C_0).^2;            %expected [Ca] background
  adj=0;
  Z=exp((adj-dC)/(2*s2))/sqrt(2*pi*s2);%[Ca] prob
  z=repmat(p,size(Z)); z(n>0)=1-p;%spike prob
  Z=Z.*z;
elseif(isempty(X1))
  n=X0(1); C=X0(2);           %unwrap previous point
  n=rand>p;                   %new spike    
  C=C-P.dt/P.tau_c*(C-P.C_0)+P.A*n+sqrt(s2)*randn; %new Ca
  Z=[n C];
else
  n0=X0(:,1); C0=X0(:,2); M=length(C0);
  n1=X1(:,1); C1=X1(:,2); N=length(C1);

  C0=C0-P.dt/P.tau_c*(C0-P.C_0);          %expected [Ca]-left
  C1=C1-P.A*n1;                           %expected [Ca]-right
  A=repmat(reshape(C1(:),N,1),[1,M,1]);
  B=repmat(reshape(C0(:),1,M),[N,1,1]);
  A=(A-B).^2;
  adj=0;
  Z=exp((adj-A)/(2*s2))/sqrt(2*pi*s2);    %[Ca]-[Ca] transition prob

  z=repmat(p,N,1); z(n1>0)=1-p;           %n1 generation prob
  Z=Z.*repmat(z(:),[1 M]);
end