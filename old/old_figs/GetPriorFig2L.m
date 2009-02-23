function GetPriorFig2L(Sim,R,S,M)

O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.K_o);
Oind    = find(~isnan(O));
xs      = [Sim.tvec(Oind(2)) 1];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

priAX   = [xs 0 1];
parAX   = [xs 0 1];
cmin1   = min(min(S(1).C(:,xind)));
cmax1   = max(max(S(1).C(:,xind)));

gray    = [0.75 0.75 0.75];

col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8;
ccol(ccol>1)=1;
ind     = Sim.K:-1:1;
Nrows   = 1;

cmin2   = min(min(S(2).C(:,xind)));
cmax2   = max(max(S(2).C(:,xind)));
cmin    = min(cmin1,cmin2);
cmax    = max(cmax1,cmax2);
Ncols   = 2;

figure(2), clf,
set(gcf, 'color', 'w');

%% backwards particle filter
l1   = .2;
w   = .35;
b   = .15;
h   = .33;
hs  = .05;
l2  = .05+l1+w;
tl  = .05;
yfs = 12;
xfs = 8;
tfs = 14;

% plot forward particles
subplot('Position',[l1 (b+h+hs) w h]) %subplot('Position',[left bottom width height])
cla, hold on
plot(Sim.tvec,(S(1).C'-cmin)/(cmax-cmin),'Color',ccol(2,:))
plot(Sim.tvec,(M(1).fCbar-cmin)/(cmax-cmin),'Color',col(2,:),'linewidth',2)
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',.5,'markersize',3)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'color',gray,'LineWidth',1)
axis(parAX)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
ylab=ylabel({'Calcium';'Particles'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs)
title('Backwards','fontsize',tfs)

subplot('Position',[l1 b w h]) %subplot('Position',[left bottom width height])
stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
BarVar1=M(1).fIbar+M(1).fIvar;
BarVar1(BarVar1>1)=1;
stem(Sim.tvec,BarVar1,'Marker','none','Color',ccol(2,:),'LineWidth',2)
stem(Sim.tvec,M(1).fIbar,'Marker','none','Color',col(2,:),'LineWidth',2)
axis(parAX)
set(gca,'YTickLabel',[])
ylab=ylabel({'Spike';'Histogram'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs)
xlab=xlabel('Time (sec)');
set(xlab,'fontsize',xfs)

%% vanilla particle filter

subplot('Position',[l2 (b+h+hs) w h]) %subplot('Position',[left bottom width height])
cla, hold on
plot(Sim.tvec,(S(2).C'-cmin)/(cmax-cmin),'Color',ccol(3,:))
plot(Sim.tvec,(M(2).fCbar-cmin)/(cmax-cmin),'Color',col(3,:),'linewidth',2)
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',0.5,'markersize',3)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'color',gray,'LineWidth',1)
axis(parAX)
set(gca,'YTickLabel',[])
set(gca,'XTickLabel',[])
title('Prior','fontsize',tfs)

subplot('Position',[l2 b w h]) %subplot('Position',[left bottom width height])
stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
BarVar2=M(2).fIbar+M(2).fIvar;
BarVar2(BarVar2>1)=1;
stem(Sim.tvec,BarVar2,'Marker','none','Color',ccol(3,:),'LineWidth',2)
stem(Sim.tvec,M(2).fIbar,'Marker','none','Color',col(3,:),'LineWidth',2)
axis(parAX)
set(gca,'YTick',[]), set(gca,'YTickLabel',[])
xlab=xlabel('Time (sec)');
set(xlab,'fontsize',xfs)

fig=figure(2);
bgr=1.2*[7 5];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior