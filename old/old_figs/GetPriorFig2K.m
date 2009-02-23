function GetPriorFig2K(Sim,R,S,M)

O       = R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
Oind    = find(~isnan(O));
xs      = [Sim.tvec(Oind(2)) Sim.tvec(Oind(end-2))];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

priAX   = [xs 0 1];
parAX   = [xs 0 2];
S1      = S{1};
M1      = M(1);
cmin1   = min(min(S1.C(:,xind)));
cmax1   = max(max(S1.C(:,xind)));
Cbar1   = sum(S1.w_f.*S1.C);
Ibar1   = sum(S1.w_f.*S1.I);
Ivar1   = sum((repmat(Ibar1,Sim.N,1)-S1.I).^2)/Sim.N;

gray    = [0.75 0.75 0.75];

col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8;
ccol(ccol>1)=1;
ind     = Sim.K:-1:1;
Nrows   = 3;
if numel(Sim.pfs)==2
    S2      = S{2};
    M2      = M(2,1);
    cmin2   = min(min(S2.C(:,xind)));
    cmax2   = max(max(S2.C(:,xind)));
    Cbar2   = 1/Sim.N*sum(S2.C);
    Ibar2   = 1/Sim.N*sum(S2.I);
    Ivar2   = sum((repmat(Ibar2,Sim.N,1)-S2.I).^2)/Sim.N;
    cmin    = min(cmin1,cmin2);
    cmax    = max(cmax1,cmax2);
    Ncols   = 2;
else
    cmin=cmin1;
    cmax=cmax1;
    Ncols=1;
end

figure(2), clf,
set(gcf, 'color', 'w');

%% backwards particle filter
l1   = .2;
w   = .35;
b   = .15;
h   = .22;
hs  = .05;
l2  = .05+l1+w;
tl  = .05;
ylfs= 12;
tfs = 14;

% plot prior
subplot('Position',[l1 b+(h+hs)*2 w h]) %subplot('Position',[left bottom width height])
filh=fill([Sim.tvec Sim.tvec(ind)],[R.p zeros(1,Sim.K)],gray);
set(filh,'edgecolor',gray)
axis(priAX)
set(gca,'XTickLabel',[]),
set(gca,'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(Oind),'TickLength',[tl .25]),
set(gca,'YTick',[])
set(gca,'box','off')
ylab=ylabel({'Probability'; 'of Spiking'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',ylfs)
title('Backwards','fontsize',tfs)
  
% plot forward particles
subplot('Position',[l1 (b+h+hs) w h]) %subplot('Position',[left bottom width height])
cla, hold on
plot(Sim.tvec,(S1.C'-cmin)/(cmax-cmin)+1,'Color',ccol(2,:))
plot(Sim.tvec,(Cbar1-cmin)/(cmax-cmin)+1,'Color',col(2,:),'linewidth',2)
plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',.5,'markersize',3)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')
stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
BarVar1=Ibar1+Ivar1;
BarVar1(BarVar1>1)=1;
stem(Sim.tvec,BarVar1,'Marker','none','Color',ccol(2,:),'LineWidth',2)
stem(Sim.tvec,Ibar1,'Marker','none','Color',col(2,:),'LineWidth',2)

axis(parAX)
set(gca,'XTickLabel',[]),
set(gca,'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(Oind),'TickLength',[tl .25])
set(gca,'YTick',[])
ylab=ylabel({'Forward';'Particles'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',ylfs) 
% plot forward backward
subplot('Position',[l1 b w h]) %subplot('Position',[left bottom width height])
cla, hold on

hfill=fill([Sim.tvec Sim.tvec(ind)],([M1.Cbar-sqrt(M1.Cvar) M1.Cbar(ind)+sqrt(M1.Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(M1.Cbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')

stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
BarVar1=M1.Ibar+M1.Ivar;
BarVar1(BarVar1>1)=1;
stem(Sim.tvec,BarVar1,'Marker','none','Color',ccol(2,:),'LineWidth',2)
stem(Sim.tvec,Ibar1,'Marker','none','Color',col(2,:),'LineWidth',2)

axis(parAX)
set(gca,'YTick',[]), set(gca,'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(Oind),'TickLength',[tl .25])
% XTickLabels=Sim.tvec(Oind)-Sim.tvec(Oind(2));
% XTickLabels={'';XTickLabels(2);'';'';'';'';XTickLabels(7);'';'';'';'';XTickLabels(12);'';'';'';''};
set(gca,'fontsize',8),
ylab=ylabel({'Inferred';'Distributions'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',ylfs)
xlab=xlabel('Time (sec)');
set(xlab,'fontsize',8)

%% vanilla particle filter
if numel(Sim.pfs)==2
    subplot('Position',[l2 b+(h+hs)*2 w h]) %subplot('Position',[left bottom width height])
    cla, hold on
    filh=fill([Sim.tvec Sim.tvec(ind)],[R.p zeros(1,Sim.K)],gray);
    set(filh,'edgecolor',gray)
    axis(priAX)
    set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
%     set(gca,'XTick',Sim.tvec(Oind),'TickLength',[tl .25]),
    set(gca,'YTick',[])
    set(gca,'box','off')
    title('Prior','fontsize',tfs)

    subplot('Position',[l2 (b+h+hs) w h]) %subplot('Position',[left bottom width height]) 
    cla, hold on
    plot(Sim.tvec,(S2.C'-cmin)/(cmax-cmin)+1,'Color',ccol(3,:))
    plot(Sim.tvec,(Cbar2-cmin)/(cmax-cmin)+1,'Color',col(3,:),'linewidth',2)
    plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',0.5,'markersize',3)
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')
    stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
    BarVar2=Ibar2+Ivar2;
    BarVar2(BarVar2>1)=1;
    stem(Sim.tvec,BarVar2,'Marker','none','Color',ccol(3,:),'LineWidth',2)
    stem(Sim.tvec,Ibar2,'Marker','none','Color',col(3,:),'LineWidth',2)

    axis(parAX)
    set(gca,'YTick',[]), set(gca,'YTickLabel',[])
%     set(gca,'XTick',Sim.tvec(Oind),'TickLength',[tl .25])
    set(gca,'XTickLabel',[]),

    subplot('Position',[l2 b w h]) %subplot('Position',[left bottom width height])
    cla, hold on
    hfill=fill([Sim.tvec Sim.tvec(ind)],([M2.Cbar-sqrt(M2.Cvar) M2.Cbar(ind)+sqrt(M2.Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(3,:));
    set(hfill,'edgecolor',ccol(3,:))
    plot(Sim.tvec,(M2.Cbar-cmin)/(cmax-cmin)+1,'linewidth',1,'color',col(3,:))
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
    BarVar1=M2.Ibar+M2.Ivar;
    BarVar1(BarVar1>1)=1;
    stem(Sim.tvec,BarVar1,'Marker','none','Color',ccol(3,:),'LineWidth',2)
    stem(Sim.tvec,Ibar1,'Marker','none','Color',col(3,:),'LineWidth',2)

    axis(parAX)
    set(gca,'YTick',[]), set(gca,'YTickLabel',[])
%     set(gca,'XTick',Sim.tvec(Oind),'TickLength',[tl .25])
%     XTickLabels=Sim.tvec(Oind)-Sim.tvec(Oind(2));
%     XTickLabels={'';XTickLabels(2);'';'';'';'';XTickLabels(7);'';'';'';'';XTickLabels(12);'';'';'';''};
     set(gca,'fontsize',8),
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    xlab=xlabel('Time (sec)');
    set(xlab,'fontsize',8)
end

fig=figure(2);
bgr=1.2*[7 5];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior