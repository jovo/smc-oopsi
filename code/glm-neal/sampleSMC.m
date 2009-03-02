function [Xr,collector,Au]=sampleSMC(P,Y,funcPY,funcPXX,funcG)
%Sequential Monte Carlo for conditioned samples in continuous HMM.
%[Xr,Tr,Au]=sampleSMC(p,Yr,@PY,@PXX,@PG) computes chain of states 
%conditioned on observations Yr in HMM specified by observation 
%probabilities PY and transition probabilities PXX. See help for sample 
%spPY and spPXX for definition of the functions setting HMM. Yr is 
%cell-array of T observations. "p" is a structure of parameters to be 
%passed onto PY and PXX, additionally defining P.iter - number of 
%iterations in SMC - and P.grid - number of points to use in stochastic 
%"grid". Au is a vector of autocorrelation values, and Tr is a 1xK
%cell-array of K sampled M dimensional states. YM(01/2009)

T=length(Y);
K=P.iter; %number of iterations of MC till convergence
N=P.grid; %size of stochastic grid to use

rej=0;
if(nargout>1) collector=cell(1,K); end %sampler history

X=cell(1,T); %form stochastic grid /pool/
for t=1:T 
  [tmp p]=funcG(P,[],t);
  X{t}=zeros(N,length(tmp)+1); X{t}(1,:)=[tmp p]; %keep rho(X,t) for ref
  for k=2:N 
    [tmp p]=funcG(P,[],t); 
    X{t}(k,:)=[tmp p];  
  end
end

for it=1:K
  [Xr xr]=sampleMC(P,X,Y,funcPY,funcPXX); %draw sample MCC
  if(it>1 && max(xr)==1) rej=rej+1; end   %record rejection  
  if(nargout>1) collector{it}=Xr; end
  
  if(it<K)
    for t=1:T %reform stochastic grid /pool/
      X{t}(1,:)=Xr{t};  %pulled chain should be in
      for k=2:N
        [tmp p]=funcG(P,[],t);
        X{t}(k,:)=[tmp p];
      end
    end
    if(mod(it,ceil(K/25))==0) fprintf('.'); end
  else
    fprintf('a.r.=%.3g\n',1-rej/(K-1));
  end
end

if(nargout>2) %obtain autocovariance measurement
  n=length(Xr{1}(:))-1;
  z=zeros(K,n*T); 
  for it=1:K 
    z1=zeros(T,n);
    for t=1:T z1(t,:)=collector{it}{t}(1:n); end
    z(it,:)=z1(:)';
  end
  z1=xcov(z(:,1),'unbiased'); Au=zeros(length(z1),n*T); Au(:,1)=z1;
  for it=2:n*T Au(:,it)=xcov(z(:,it),'unbiased'); end
end