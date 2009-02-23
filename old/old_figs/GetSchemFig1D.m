function GetSchemFig1D(Sim,R,S,B,M)

M       = M(1,1);
S       = S{1};
O       = R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
xs      = [Sim.tvec(35) Sim.tvec(135)];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

pmax = max(R.p(xind));
pAX  = [xs 0 pmax];

spikemax= max(1,max(R.I(xind)));
spikeAX = [xs 0 spikemax];

cmin=min(min(min(R.C(xind)),min(M.Cbar(xind)-sqrt(M.Cvar(xind)))),min(O(xind)));
cmax=max(max(max(R.C(xind)),max(M.Cbar(xind)+sqrt(M.Cvar(xind)))),max(O(xind)));

gray=[0.75 0.75 0.75];

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;
ind=Sim.K:-1:1;

figure(1), clf, Nsubs=5;

%% filtered stimulus
i=1; subplot(Nsubs,1,i), cla, hold on
stem(Sim.tvec,B.kx,'marker','none','color','k','LineWidth',8)
plot(.1,10)
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
%set(gca,'XTick',[]), 
set(gca,'YTick',[])
axis([xs min(B.kx(xind)) max(B.kx(xind))])
ylab=ylabel({'External';'Stimulus'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% hidden states
i=i+1; subplot(Nsubs,1,i), cla, hold on
% plot(Sim.tvec,-1*R.I+1.95,'+','Color',gray,'LineWidth',2)
stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray,'LineWidth',2)
% plot(Sim.tvec,R.p/pmax,'-','Color',gray,'LineWidth',1)
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
%set(gca,'XTick',[]), 
set(gca,'YTick',[])
axis(spikeAX)
ylab=ylabel({'[Ca^2^+] and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% observation state
i=i+1; subplot(Nsubs,1,i), cla, hold on
% stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
% plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray)
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',2)
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
%set(gca,'XTick',[]), 
set(gca,'YTick',[])
axis(spikeAX)
ylab=ylabel({'Observations'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% inferred calcium
i=i+1; subplot(Nsubs,1,i), cla, hold on, 

h=fill([Sim.tvec Sim.tvec(ind)],([M.Cbar-sqrt(M.Cvar) M.Cbar(ind)+sqrt(M.Cvar(ind))]-cmin)/(cmax-cmin),ccol(1,:));
set(h,'edgecolor',ccol(1,:))
plot(Sim.tvec,(M.Cbar-cmin)/(cmax-cmin),'linewidth',2,'color',col(1,:))
% plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',2)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'color',gray,'LineWidth',2)

set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
%set(gca,'XTick',[]), 
set(gca,'YTick',[])
axis(spikeAX)
ylab=ylabel({'Calcium';'Distribution'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% inferred spikes
i=i+1; subplot(Nsubs,1,i), cla, hold on, 

stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
stem(Sim.tvec,M.Ibar+M.Ivar,'Marker','none','Color',ccol(1,:),'LineWidth',2)
stem(Sim.tvec,M.Ibar,'Marker','none','Color',col(1,:),'LineWidth',2)

set(gca,'YTickLabel',[]), %set(gca,'XTickLabel',[])
set(gca,'YTick',[]), %set(gca,'XTick',[])
axis(spikeAX)
ylab=ylabel({'Spike';'Histogram'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlabel('Time (ms)')

fig=figure(1);
bgr=0.5*[7 7];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
%print -depsc C:\D\Research\liam\SMC_EM_GLM\schem;