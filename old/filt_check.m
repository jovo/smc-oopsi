dt=.01;
t=0:dt:5;
tau=.3;

kern=exp(-t/tau);
s=zeros(1,length(t));
s(100)=1;

ks1=conv(kern,s);
ks2=filter(1,[1 -1/tau],s);

plot(ks1),
hold on
plot(ks2,'r')
ylabel(max(ks2))


for i=1:10000;
ass1=filter(1,[1 -(1-.01/tau)],s);
ass2=conv(kern,s);
end
