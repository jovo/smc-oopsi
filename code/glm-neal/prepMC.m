function [Yr,Xr]=prepMC(P,X,funcPY,funcPXX)
%Sequence of states and observation for HMM.
%[Yr,Xr]=prepMC(p,X,@PY,@PXX) computes sequence of observations Yr and 
%either discrete or continuous states Xr for HMM specified by observation 
%probabilities PY and transition probabilities PXX. See help for sample 
%spPY and spPXX for definition of the functions setting HMM. If X is a 
%number, prepMC computes a sequence of continuous states Xr(t), as returned 
%by PXX, of length X=T. If X is a cell-array of NxM matrices of length T, 
%prepMC computes a sequence of states such that at each time Xr is one of 
%N M-dimensional states in X{t}. Both Yr and Xr are 1xT cell-arrays 
%containing instances of multidimensional observations and states at each 
%time point t. "p" is a structure of parameters to be passed onto PY and 
%PXX. YM(01/2009)

T=length(X);
Yr=cell(1,T);
Xr=cell(1,T);

if(length(X)==1) %if continuous
  T=X;
  %INITIALIZE STATES P(x_1|Y_1)=P(Y_1|x_1)/Z
  Xr{1}=funcPXX(P,[]); Yr{1}=funcPY(P,[],Xr{1},1);

  %FORWARD PASS  
  for t=2:T Xr{t}=funcPXX(P,[],Xr{t-1},t); Yr{t}=funcPY(P,[],Xr{t},t); end
else
  %INITIALIZE STATES P(x_1|Y_1)=P(Y_1|x_1)/Z
  p=funcPXX(P,X{1}); p=cumsum(p)/sum(p); %normalization
  r=rand; i=find(p>r,1); 
  Xr{1}=X{1}(i,:); Yr{1}=funcPY(P,[],Xr{1},1);

  %FORWARD PASS
  for t=2:T
    p=funcPXX(P,X{t},Xr{t-1},t-1); p=cumsum(p)/sum(p); %normalization
    r=rand; i=find(p>r,1); 
    Xr{t}=X{t}(i,:); Yr{t}=funcPY(P,[],Xr{t},t);
  end
end