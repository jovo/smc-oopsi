clc,
spikemat=zeros(2^Sim.freq,Sim.freq);
B.beta=2;
s=5;
for tt=s:-1:1, 
%     spikemat(1:2^(s-tt+1),tt)=repmat([0 2],1,2^(s-tt))';
    spikemat(1:2^(s-tt+1),tt)=[zeros(1,2^(s-tt)) B.beta*ones(1,2^(s-tt))]';
end
spikemat