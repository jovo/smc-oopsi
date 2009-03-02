function [Z p]=PG_GOOPSI(P,X,t)
%Grid function rho(X,t) [for GOOPSI:Vogelstein-Paninski setup]
%Z=PG(P,[],t) returns a sample of X from rho(X,t).
%Z=PG(p,X,t) computes rho(X,t). If X is a NxM matrix of N M-dimensional
%states, Z is Nx1 vector of P_i=rho(X_i,t).
%[Z,P]=PG(p,[],t) returns a sample of X from rho(X,t) and P=rho(X,t).
%"p" is a structure of parameters. YM(01/2009)
%
%REQUIRED PARAMETERS
% P.nf           1xT cell array of particles.
% P.wf           MxT cell array of particle weights.
% P.vf           2xT array of kernel variances.

if(isempty(X))
  %THIS IS THE ONLY PART THAT SHOULD RUN WITH SAMPLER  
  ph=P.wf(:,t);             %weights of different gaussians
  phc=cumsum(ph)/sum(ph);   %normalization
  r=rand; i=find(r<phc,1);  %draw the base state
  
  n=P.nf{t}(i,1); C=P.nf{t}(i,2);%unwrap the base state
  
  n=(rand>P.vf(1,t))==(n>0);%randomize spike-state, trick here
  C=C+P.vf(2,t)*randn;      %randomize Ca-state
  Z=[n C];

  if(nargout==1) p=NaN; return; end
                            %evaluate rho(X,t)=P(n)*P(C)
  p1=sum(ph(P.nf{t}(:,1)>0));%prob to pull a particle with n==1
  if(n>0)                   %if n==1, P(n)=P(base=1)*P{stay}+P(base=0)*P{switch}
    p=(1-P.vf(1,t))*p1+P.vf(1,t)*(1-p1);
  else                      %if n==0, P(n)=P(base=1)*P{switch}+P(base=0)*P{stay}
    p=P.vf(1,t)*p1+(1-P.vf(1,t))*(1-p1);
  end
  dC=(C-P.nf{t}(:,2)).^2;   %prob to pull C from diff particles
  s2=P.vf(2,t)^2;
  p1=exp(-dC/(2*s2));
  p=p*sum(ph.*p1);
else
  Z=zeros(size(X,1),1);
  for k=1:length(Z)
    n=X(k,1); C=X(k,2);     %unwrap state
    p1=sum(ph(P.nf{t}(:,1)>0));%prob to pull particle with n=1
    if(n>0)                 %if n==1, P(n)=P(base=1)*P{stay}+P(base=0)*P{switch}
      Z(k)=(1-P.vf(1,t))*p1+P.vf(1,t)*(1-p1);
    else                    %if n==0, P(n)=P(base=1)*P{switch}+P(base=0)*P{stay}
      Z(k)=P.vf(1,t)*p1+(1-P.vf(1,t))*(1-p1);
    end
    dC=(C-P.nf{t}(:,2)).^2; %prob to pull C from diff particles
    s=P.v(2,t)^2;
    p1=exp(-dC/(2*s));
    Z(k)=Z(k)*sum(ph.*p1);
  end
  p=NaN;
end


function C = AHill_v1(P,F)
% generalized hill model, inverse
C = (F./(1-F)).^(1/P.n)*P.k_d; C(C<0)  = 0;
