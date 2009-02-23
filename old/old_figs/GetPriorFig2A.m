function GetPriorFig2A(Sim,R,S,B,M)

O       = R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
S1=S{2}; S2=S{6};
M1=M(2,1); M2=M(6,1);
xs      = [1:find(Sim.tvec>0.5,1)];
AX      = [min(Sim.tvec(xs)) max(Sim.tvec(xs)) -1 1];
cmin1   = min(min(S1.C(:,xs)));
cmin2   = min(min(S2.C(:,xs)));
cmax1   = max(max(S1.C(:,xs)));
cmax2   = max(max(S2.C(:,xs)));
gray    = [0.75 0.75 0.75];

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;
ind=Sim.K:-1:1;

figure, clf, 

%%
subplot(2,2,1), cla, hold on
plot(Sim.tvec,B.kx,'k','LineWidth',2)
set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',[]), set(gca,'YTick',[])
axis([min(Sim.tvec(xs)) max(Sim.tvec(xs)) min(B.kx(xs)) max(B.kx(xs))])

subplot(2,2,2), cla, hold on
plot(Sim.tvec,B.kx,'k','LineWidth',2)
set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',[]), set(gca,'YTick',[])
axis([min(Sim.tvec(xs)) max(Sim.tvec(xs)) min(B.kx(xs)) max(B.kx(xs))])

%%
subplot(2,2,3), cla, hold on
plot(Sim.tvec,(S1.C'-cmin1)/(cmax1-cmin1),'Color','b')
plot(Sim.tvec,(O-cmin1)/(cmax1-cmin1),'ok','LineWidth',2)
plot(Sim.tvec,(R.C-cmin1)/(cmax1-cmin1),'color',gray,'LineWidth',2)
h=fill([Sim.tvec Sim.tvec(ind)],[M1.Ibar-M1.Ivar M1.Ibar(ind)+M1.Ivar(ind)]-1,ccol(1,:));
set(h,'edgecolor',ccol(1,:))
plot(Sim.tvec,M1.Ibar-1,'linewidth',2,'color',col(1,:))
stem(Sim.tvec,-1*R.I,'Marker','none','Color',gray,'LineWidth',2)
axis(AX)
set(gca,'YTickLabel',[])

%%
subplot(2,2,4), cla, hold on
plot(Sim.tvec,(S2.C'-cmin2)/(cmax2-cmin2),'Color','b')
plot(Sim.tvec,(O-cmin2)/(cmax2-cmin2),'ok','LineWidth',2)
plot(Sim.tvec,(R.C-cmin2)/(cmax2-cmin2),'color',gray,'LineWidth',2)
h=fill([Sim.tvec Sim.tvec(ind)],[M2.Ibar-M2.Ivar M2.Ibar(ind)+M2.Ivar(ind)]-1,ccol(2,:));
set(h,'edgecolor',ccol(2,:))
plot(Sim.tvec,M2.Ibar-1,'linewidth',2,'color',col(2,:))
stem(Sim.tvec,-1*R.I,'Marker','none','Color',gray,'LineWidth',2)
axis(AX)
set(gca,'YTickLabel',[])
