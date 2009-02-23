function GetPriorFig2J(Sim,R,S,M)

O       = R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
Oind    = find(~isnan(O));
xs      = [Sim.tvec(Oind(2)) Sim.tvec(Oind(end-1))];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

priAX   = [xs 0 1];
parAX   = [xs 0 2];
S1      = S{1};
M1      = M(1,1);
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
% plot prior
subplot(Nrows,Ncols,1)
h=fill([Sim.tvec Sim.tvec(ind)],[R.p zeros(1,Sim.K)],gray);
set(h,'edgecolor',gray)
axis(priAX)
set(gca,'XTickLabel',[]),
set(gca,'YTickLabel',[])
set(gca,'XTick',Sim.tvec(Oind)),
set(gca,'YTick',[])
ylab=ylabel({'Probability'; 'of Spiking'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',8)
title('Backwards')
 
% plot forward particles
subplot(Nrows,Ncols,3), 
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
set(gca,'XTick',Sim.tvec(Oind)),
set(gca,'YTick',[])
ylab=ylabel({'Forward';'Particles'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',8)

% plot forward backward
subplot(Nrows,Ncols,5), cla, hold on

h=fill([Sim.tvec Sim.tvec(ind)],([M1.Cbar-sqrt(M1.Cvar) M1.Cbar(ind)+sqrt(M1.Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
set(h,'edgecolor',ccol(2,:))
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
set(gca,'XTick',Sim.tvec(Oind)),
XTickLabels=Sim.tvec(Oind)-Sim.tvec(Oind(2));
XTickLabels={'';XTickLabels(2);'';'';'';'';XTickLabels(7);'';'';'';'';XTickLabels(12);'';'';'';''};
set(gca,'XTickLabel',XTickLabels,'fontsize',8),
ylab=ylabel({'Inferred';'Distributions'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',8)
xlab=xlabel('Time (sec)');
set(xlab,'fontsize',8)

% vanilla particle filter
if numel(Sim.pfs)==2
    subplot(Nrows,Ncols,2), cla, hold on
    h=fill([Sim.tvec Sim.tvec(ind)],[R.p zeros(1,Sim.K)],gray);
    set(h,'edgecolor',gray)
    axis(priAX)
    set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
    set(gca,'XTick',[]), set(gca,'YTick',[])
    title('Prior')

    subplot(Nrows,Ncols,4), cla, hold on
    plot(Sim.tvec,(S2.C'-cmin)/(cmax-cmin)+1,'Color',ccol(3,:))
    plot(Sim.tvec,(Cbar2-cmin)/(cmax-cmin)+1,'Color',col(3,:),'linewidth',2)
    plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',0.5,'markersize',3)
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')
    stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
    BarVar2=Ibar2+Ivar2;
    BarVar2(BarVar2>1)=1;
    stem(Sim.tvec,BarVar2,'Marker','none','Color',ccol(3,:),'LineWidth',2)
    stem(Sim.tvec,Ibar2,'Marker','none','Color',col(3,:),'LineWidth',2)

    axis(parAX)
    set(gca,'YTick',[]), set(gca,'YTickLabel',[])
    set(gca,'XTick',Sim.tvec(Oind)),
    set(gca,'XTickLabel',[]),

    subplot(Nrows,Ncols,6), cla, hold on
    h=fill([Sim.tvec Sim.tvec(ind)],([M2.Cbar-sqrt(M2.Cvar) M2.Cbar(ind)+sqrt(M2.Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(3,:));
    set(h,'edgecolor',ccol(3,:))
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
    set(gca,'XTick',Sim.tvec(Oind)),
    XTickLabels=Sim.tvec(Oind)-Sim.tvec(Oind(2));
    XTickLabels={'';XTickLabels(2);'';'';'';'';XTickLabels(7);'';'';'';'';XTickLabels(12);'';'';'';''};
    set(gca,'XTickLabel',XTickLabels,'fontsize',8),
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    xlab=xlabel('Time (sec)');
    set(xlab,'fontsize',8)
end

fig=figure(2);
bgr=0.5*[7 5];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior;