function GetSchemFig1Cb(Sim,R,S,B,M)

M=M(2,1);
S=S{2};

O=R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
for k=1/Sim.dt:Sim.K-1/Sim.dt
    sumn(k)=sum(R.I(k:k+1/Sim.dt));
end
[maxn ind]=max(sumn);
int     = 0.3;
xs      = [Sim.tvec(Sim.frac*2) 0.5-Sim.dtSample];%[ind*Sim.dt-int ind*Sim.dt+int];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

pmax = max(R.p(xind));
pAX  = [xs 0 pmax];

spikemax= max(1,max(R.I(xind)));
spikeAX = [xs 0 spikemax];

cmin=min(min(min(R.C(xind)),min(M.Cbar(xind))),min(O(xind)));
cmax=max(max(max(R.C(xind)),max(M.Cbar(xind))),max(O(xind)));

gray=[0.75 0.75 0.75];

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;
ind=Sim.K:-1:1;

figure(1), Nsubs=5;

%% filtered stimulus
i=1; subplot(Nsubs,1,i), cla, hold on
plot(Sim.tvec,B.kx,'k','LineWidth',2)
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
ylab=ylabel({'Hidden';'States'});
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
ylab=ylabel({'Observation';'State'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% inferred calcium
i=i+1; subplot(Nsubs,1,i), cla, hold on, 

plot(Sim.tvec,(S.C'-cmin)/(cmax-cmin),'Color',ccol(1,:))
plot(Sim.tvec,(M.Cbar-cmin)/(cmax-cmin),'Color',col(1,:),'linewidth',2)
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',2)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'color',gray,'LineWidth',2)

set(gca,'XTickLabel',[]), 
set(gca,'YTickLabel',[])
%set(gca,'XTick',[]), 
set(gca,'YTick',[])
axis(spikeAX)
ylab=ylabel({'Calcium';'Particles'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% inferred spikes
i=i+1; subplot(Nsubs,1,i), cla, hold on, 

h=fill([Sim.tvec Sim.tvec(ind)],[M.Ibar-M.Ivar M.Ibar(ind)+M.Ivar(ind)],ccol(1,:));
set(h,'edgecolor',ccol(1,:))
plot(Sim.tvec,M.Ibar,'linewidth',1,'color',col(1,:))

% plot(Sim.tvec,-1*R.I+1.95,'+','Color',gray,'linewidth',2)
stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
% plot(Sim.tvec,(O-cmin)/max(O(xind)),'o','LineWidth',2,'Color','k');
set(gca,'YTickLabel',[]), %set(gca,'XTickLabel',[])
set(gca,'YTick',[]), %set(gca,'XTick',[])
axis(spikeAX)
ylab=ylabel({'Spike';'Distribution'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlabel('Time (ms)')

fig=figure(1);
bgr=0.5*[7 7];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\schem;
