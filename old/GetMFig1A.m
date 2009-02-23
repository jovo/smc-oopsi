function GetMFig1A(Sim,R,S,M,E)

O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.K_o);
Oind    = find(~isnan(O));
xs      = [Sim.tvec(Oind(2)) Sim.tvec(Oind(end-1))];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

spikeAX = [xs 0 1];
hidAX   = [xs 0 2];

cmin=min(min(min(R.C(xind)),min(M.bCbar(xind)-sqrt(M.bCvar(xind)))),min(O(xind)));
cmax=max(max(max(R.C(xind)),max(M.bCbar(xind)+sqrt(M.bCvar(xind)))),max(O(xind)));

gray=[0.75 0.75 0.75];

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;
ind=Sim.K:-1:1;

figure(1), clf, Nsubs=2;
set(gcf, 'color', 'w');

%% linear kernel
i=1; subplot(Nsubs,1,i), cla, hold on
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
set(gca,'XTick',Sim.tvec(Oind)), 
set(gca,'YTick',[])
axis([xs min(Sim.x(xind)) max(Sim.x(xind))+1.])
ylab=ylabel({'External';'Stimulus'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% other parameters
i=i+1; subplot(Nsubs,1,i), cla, hold on
set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
set(gca,'XTick',Sim.tvec(Oind)), 
set(gca,'YTick',[])
axis(hidAX)
ylab=ylabel({'Calcium';'and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',gray)

fig=figure(1);
wh=[6 3];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\schem;