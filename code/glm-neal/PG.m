function Z=PG(P,X,t)
%Z=PG(p,X,t)
%Grid function rho(X,t). If X is empty, returns a sample of X from 
%rho(X,t). In case X is a vector, produces a vector of P_i=rho(X_i,t). "p"
%is a structure of parameters.
if(isempty(X))
  Z=randn*sqrt(t);
else  
  Z=exp(-X.^2/2/t)/sqrt(2*pi*t);
end
