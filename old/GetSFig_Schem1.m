function GetSchemFig1G(Sim,R,S,M)

O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.K_o);
Oind    = find(~isnan(O));
xs      = [Sim.tvec(Oind(2)) Sim.tvec(Oind(end-1))];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

pmax = max(R.p(xind));
pAX  = [xs 0 pmax];

spikeAX = [xs 0 1];
hidAX   = [xs 0 2];

cmin=min(min(min(R.C(xind)),min(M.bCbar(xind)-sqrt(M.bCvar(xind)))),min(O(xind)));
cmax=max(max(max(R.C(xind)),max(M.bCbar(xind)+sqrt(M.bCvar(xind)))),max(O(xind)));

gray=[0.75 0.75 0.75];

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;
ind=Sim.K:-1:1;

figure(1), clf, Nsubs=4;
set(gcf, 'color', 'w');

%% external stimulus
i=1; subplot(Nsubs,1,i), cla, hold on
plot(Sim.tvec,Sim.x,'k','LineWidth',2)
plot(.1,10)
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
set(gca,'XTick',Sim.tvec(Oind)), 
set(gca,'YTick',[])
axis([xs min(Sim.x(xind)) max(Sim.x(xind))+1.])
ylab=ylabel({'External';'Stimulus'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% hidden states
i=i+1; subplot(Nsubs,1,i), cla, hold on
stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'Color',gray,'LineWidth',2)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
set(gca,'XTick',Sim.tvec(Oind)), 
set(gca,'YTick',[])
axis(hidAX)
ylab=ylabel({'Calcium';'and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',gray)

%% observation state
i=i+1; subplot(Nsubs,1,i), cla, hold on
% stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
% plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray)
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',1,'markersize',4)
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
set(gca,'XTick',Sim.tvec(Oind)), 
set(gca,'YTick',[])
axis(spikeAX)
ylab=ylabel({'Observations'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% inferred calcium
i=i+1; subplot(Nsubs,1,i), cla, hold on, 

h=fill([Sim.tvec Sim.tvec(ind)],([M.bCbar-sqrt(M.bCvar) M.bCbar(ind)+sqrt(M.bCvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
set(h,'edgecolor',ccol(2,:))
plot(Sim.tvec,(M.bCbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')

stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
BarVar=M.bIbar+M.bIvar;
BarVar(BarVar>1)=1;
stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',2)
stem(Sim.tvec,M.bIbar,'Marker','none','Color',col(2,:),'LineWidth',2)

set(gca,'YTick',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',Sim.tvec(Oind)), 
set(gca,'XTickLabel',{'';'0';'';'0.1';'';'0.2';'';'0.3';'';'0.4';'';'0.5';''})%Sim.tvec(Oind)-Sim.tvec(Oind(2))), 
axis(hidAX)
ylab=ylabel({'Calcium';'and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(2,:))
xlab=xlabel('Time (sec)');

fig=figure(1);
wh=[6 3];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\schem;