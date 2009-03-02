function [Xr,xr]=sampleMC(P,X,Y,funcPY,funcPXX)
%Forward-backward procedure for conditioned samples in discrete HMM.
%[Xr,xr]=sampleMC(p,X,Yr,@PY,@PXX) computes chain of states conditioned 
%on observations Yr in HMM specified by observation probabilities PY and 
%transition probabilities PXX. See help for sample spPY and spPXX for 
%definition of the functions setting HMM. X is a cell-array of length T 
%of NxM matrices specifying at each time t the N M-dimensional states in 
%discrete states-space. Yr is cell-array of T observations. Xr is the chain
%of states, xr is the chain of integer references identifying Xr in X. 
%"p" is a structure of parameters to be passed onto PY and PXX. YM(01/2009)

T=length(Y);
N=size(X{1},1); for t=2:T N=max(N,size(X{t},1)); end
p=zeros(N,T);
pf=zeros(N,T);

%INITIALIZE STATES P(x_1|Y_1)=p(y_1|x_1)p(x_1)/Z
nh=size(X{1},1);
p(1:nh,1)=funcPY(P,Y{1},X{1},1).*funcPXX(P,X{1}); %P(x_1|Y_1)
p(1:nh,1)=p(1:nh,1)/sum(p(1:nh,1));   %normalization

%FORWARD PASS - obtain P(x_T|Y_T)
for t=2:T
  %P(x_t|Y_t)=P(y_t|x_t)*int\dx_{t-1}P(x_t|x_{t-1})P(x_{t-1}|Y_{t-1})/Z
  nh=size(X{t},1); nhprev=size(X{t-1},1);  
  p(1:nh,t)=funcPY(P,Y{t},X{t},t);
  pf(1:nh,t)=sum(funcPXX(P,X{t},X{t-1},t-1).*repmat(p(1:nhprev,t-1)',nh,1),2);  
  
  p(1:nh,t)=p(1:nh,t).*pf(1:nh,t); 
  if(sum(p(1:nh,t))==0) error('Degenerate probability field!'); end  
  p(1:nh,t)=p(1:nh,t)/sum(p(1:nh,t)); %normalization
end

%BACKWARD PASS
ph=cumsum(p(:,T))/sum(p(:,T)); 
r=rand; i=find(r<ph,1);               %produce s(x_T|Y_T)
Xr=cell(1,T); Xr{T}=X{T}(i,:);        %sample state vector
xr=zeros(1,T); xr(T)=i;               %sample state id
for t=T-1:-1:1
  %s(x_t|Y_T,x_{t+1},...)=P(x_t|Y_t)*P(x_{t+1}|x_t)/P(x_{t+1}|Y_t)
  nh=size(X{t},1); nhnext=size(X{t+1},1);    
  if(pf(i,t+1)==0) error('Error: P(x_{t+1}|Y_t)==0 in backward pass!\n'); end
  ph=funcPXX(P,Xr{t+1},X{t},t).*p(1:nh,t)'/pf(xr(t+1),t+1)'; 
  if(sum(ph)==0) error('Degenerate probability field!'); end
  ph=cumsum(ph)/sum(ph); 
  r=rand; i=find(r<ph,1);             %sample s(x_t|Y_T,x_{t+1},...)  
  xr(t)=i; Xr{t}=X{t}(i,:);
end
