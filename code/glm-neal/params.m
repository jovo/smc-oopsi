%REQUIRED PARAMETERS
P=[];

P.alpha=1; %fluorescence amplitude
P.beta=0.1; %fluorescence background
P.gamma=1e-3; %fluorescence conversion factor/noise weight
P.delta=0.01;  %background noise
P.n=1; %fluorescence saturation exponent
P.k_d=5; %fluorescence saturation offset

P.dt=0.01;  %time bin
P.f=10; %spike generation frequency
P.tau=0.2;  %lag of Ca reporter
P.C0=0; %Ca background
P.A=1; %Ca spike height
P.eps=sqrt(0.1*P.dt); %Ca noise height

P.iter=100; %SMC iterations in Neal's method
P.grid=25; %points in stochastic grid
P.invs=0.1;  %min points dispersion

P.mismatch_penalty=0.05; %time-error to call a mismatch for bipartite
P.expected_slack=0.1; %fraction of expected slack matches for bipartite