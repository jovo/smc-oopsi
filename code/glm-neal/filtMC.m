function p=filtMC(P,X,Y,funcPY,funcPXX)
%Forward-backward procedure for deriving posterior probabilities in HMM.
%P=filtMC(p,X,Yr,PY,PXX) computes posterior probabilities of states 
%P(X(t)|Y) in discrete-state HMM specified by specified by observation 
%probabilities PY, transition probabilities PXX and state space X. See help 
%for sample spPY and spPXX for definition of the functions setting HMM. X 
%is a cell-array of length T of NxM matrices specifying at each time t the 
%N M-dimensional states in discrete states-space. Yr is cell-array of T 
%observations. "p" is a structure of parameters to be passed onto PY and
%PXX. YM(01/2009)
T=length(Y);
N=size(X{1},1); for t=2:T N=max(N,size(X{t},1)); end
p=zeros(N,T);
pf=zeros(N,T);

%INITIALIZE STATES P(x_1|Y_1)=P(Y_1|x_1)/Z
nh=size(X{1},1);
p(1:nh,1)=funcPY(P,Y{1},X{1},1); p(1:nh,1)=p(1:nh,1)/sum(p(1:nh,1)); %normalization

%FORWARD PASS
for t=2:T
  %P(x_t|Y_t)=P(x_t|y_t)*int\dx_{t-1}P(x_t|x_{t-1})P(x_{t-1}|Y_{t-1})/Z
  nh=size(X{t},1); nhprev=size(X{t-1},1);
  p(1:nh,t)=funcPY(P,Y{t},X{t},t);
  pf(1:nh,t)=sum(funcPXX(P,X{t},X{t-1},t-1).*repmat(p(1:nhprev,t-1)',nh,1),2);  
  
  p(1:nh,t)=p(1:nh,t).*pf(1:nh,t); 
  if(sum(p(1:nh,t))==0) error('Degenerate prob field!'); end  
  p(1:nh,t)=p(1:nh,t)/sum(p(1:nh,t)); %normalization
end

%BACKWARD PASS
for t=T-1:-1:1
  %P(x_t|Y_T)=P(x_t|Y_t)*int\dx_{t+1}P(x_{t+1}|x_t)P(x_{t+1}|Y_T)/P(x_{t+1}|Y_t)
  nh=size(X{t},1); nhnext=size(X{t+1},1);  
  a=pf(1:nhnext,t+1); a1=p(1:nhnext,t+1); a=a1(a>0)./a(a>0); %careful div avoid zeros
  h=sum(funcPXX(P,X{t+1},X{t},t).*repmat(a,1,nh),1)';
  p(1:nh,t)=p(1:nh,t).*h;
end
