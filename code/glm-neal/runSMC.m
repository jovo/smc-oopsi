%RUN Neal's algorithm
P=[];
T=100; %spike train length

P.alpha=1; %fluorescence amplitude
P.beta=0.1; %fluorescence background
P.gamma=1e-4; %fluorescence conversion factor/noise weight
P.delta=0.01;  %background noise
P.n=1; %fluorescence saturation exponent
P.k_d=5; %fluorescence saturation offset

P.dt=0.01;  %time bin
P.f=30; %spike generation frequency
P.tau=0.2;  %lag of Ca reporter
P.C0=0; %Ca background
P.A=1; %Ca spike height
P.eps=sqrt(0.1*P.dt); %Ca noise height

P.iter=30; %SMC iterations in Neal's method
P.grid=25; %points in stochastic grid
P.invs=0.1;  %min points dispersion

P.mismatch_penalty=5; %time-error to call a mismatch for bipartite
P.expected_slack=0.1; %fraction of expected slack matches for bipartite

for k=1:1; %trial index
  [Yr,Xr]=prepMC(P,T,@spPY,@spPXX);  %obtain sample
  Xr{end}(1)=0; %remove last moment spike
  P.F=Yr;

  n=zeros(1,T);
  C=zeros(1,T);
  F=zeros(1,T);
  for i=1:T C(i)=Xr{i}(2); n(i)=Xr{i}(1); F(i)=Yr{i}; end;
  figure(k),clf,plot(C), hold on, stem(n,'r'),plot(F*P.k_d,'g','LineWidth',2)

  [Er Au]=sampleSMC(P,Yr,@spPY,@spPXX,@spPG);
  n1=zeros(1,length(Yr));
  C1=zeros(1,length(Yr));
  for i=1:length(Yr) C1(i)=Er{i}(2); n1(i)=Er{i}(1); end;
  figure(k),stem(0.9*n1,'c'),plot(C1,'m'),pause(0.1);

  tgt{k}=n;
  est{k}=n1;
end

figure,plot(sum(Au(:,1:100),2)),title('Neal aucov spikes')
figure,plot(sum(Au(:,101:200),2)),title('Neal aucov Ca')
rep=errorsSMC(P,est,tgt);
rep