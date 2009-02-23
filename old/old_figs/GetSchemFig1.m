function GetSchemFig1(Sim,R,B,M)

O=R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
for k=1/Sim.dt:Sim.K-1/Sim.dt
    sumn(k)=sum(R.I(k:k+1/Sim.dt));
end
[maxn ind]=max(sumn);
xs = [ind*Sim.dt-0.75 ind*Sim.dt+0.25];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind = [xmin:xmax];

pmax = max(R.p(xind));
pAX  = [xs 0 pmax];

spikemax= max(1,max(R.I(xind)));
spikeAX = [xs 0 spikemax];

cmin=min(min(min(R.C(xind)),min(M.Cbar(xind))),min(O(xind)));
cmax=max(max(max(R.C(xind)),max(M.Cbar(xind))),max(O(xind)));

gray=[0.75 0.75 0.75];

figure(1)
Nsubs=5;

%% filtered stimulus
i=1; subplot(Nsubs,1,i)
plot(Sim.tvec,B.kx,'k','LineWidth',2)
set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',[]), set(gca,'YTick',[])
axis([xs min(B.kx(xind)) max(B.kx(xind))])

%% hidden states
i=i+1; subplot(Nsubs,1,i), hold off, 
plot(Sim.tvec,-1*R.I+1.95,'+','Color','k','LineWidth',2)
hold on
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'k','LineWidth',2)
plot(Sim.tvec,R.p/pmax,'Color',gray,'LineWidth',1)
set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',[]), set(gca,'YTick',[])
axis(spikeAX)

%% observation state
i=i+1; subplot(Nsubs,1,i), hold off, 
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray)
hold on
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',2)
set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',[]), set(gca,'YTick',[])
axis(spikeAX)

%% inferred calcium
i=i+1; subplot(Nsubs,1,i), hold off, 
plot(Sim.tvec,(M.Cbar-cmin)/(cmax-cmin),'b','LineWidth',2); 
hold on
plot(Sim.tvec,(M.Cbar+sqrt(M.Cvar)/2-cmin)/(cmax-cmin),'b');
plot(Sim.tvec,(M.Cbar-sqrt(M.Cvar)/2-cmin)/(cmax-cmin),'b');

plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'LineWidth',1,'Color','k')
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'o','LineWidth',1,'Color','k');

set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',[]), set(gca,'YTick',[])
axis(spikeAX)

%% inferred spikes
i=i+1; subplot(Nsubs,1,i), hold off, 
plot(Sim.tvec,M.Ibar,'b','LineWidth',2)
hold on
plot(Sim.tvec,M.Ibar+sqrt(M.Ivar),'b','LineWidth',1);
plot(Sim.tvec,-1*R.I+1.95,'+','Color','k')%gray)
plot(Sim.tvec,(O-cmin)/max(O(xind)),'o','LineWidth',1,'Color','k');
set(gca,'YTickLabel',[]), %set(gca,'XTickLabel',[])
set(gca,'YTick',[]), %set(gca,'XTick',[])
axis(spikeAX)