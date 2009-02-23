function GetPriorFig2H(Sim,R,S)

O       = R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
xs      = [35:135];
AX      = [min(Sim.tvec(xs)) max(Sim.tvec(xs)) 0 2];
S1      = S{1};
cmin1   = min(min(S1.C(:,xs)));
cmax1   = max(max(S1.C(:,xs)));
Cbar1   = sum(S1.w_f.*S1.C);
Ibar1   = sum(S1.w_f.*S1.I);
Ivar1   = sum((repmat(Ibar1,Sim.N,1)-S1.I).^2)/Sim.N;


gray    = [0.75 0.75 0.75];

col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8;
ccol(ccol>1)=1;
ind     = Sim.K:-1:1;
Nrows   = 2;
if numel(Sim.pfs)==2
    S2      = S{2};
    cmin2   = min(min(S2.C(:,xs)));
    cmax2   = max(max(S2.C(:,xs)));
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

%%
subplot(Nrows,Ncols,1), cla, hold on
h=fill([Sim.tvec Sim.tvec(ind)],[R.p zeros(1,Sim.K)],gray);
set(h,'edgecolor',gray)
axis([min(Sim.tvec(xs)) max(Sim.tvec(xs)) 0 1])
set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
set(gca,'XTick',[]), set(gca,'YTick',[])
ylab=ylabel('P(spike)');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
title('Backwards Sampler')

%%
subplot(Nrows,Ncols,1+Ncols), cla, hold on
plot(Sim.tvec,(S1.C'-cmin)/(cmax-cmin)+1,'Color',ccol(2,:))
plot(Sim.tvec,(Cbar1-cmin)/(cmax-cmin)+1,'Color',col(2,:),'linewidth',2)
plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',1,'markersize',4)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')
stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
BarVar1=Ibar1+Ivar1;
BarVar1(BarVar1>1)=1;
stem(Sim.tvec,BarVar1,'Marker','none','Color',ccol(2,:),'LineWidth',2)
stem(Sim.tvec,Ibar1,'Marker','none','Color',col(2,:),'LineWidth',2)

axis(AX)
%set(gca,'XTickLabel',[]),
set(gca,'YTickLabel',[])
%set(gca,'XTick',[]),
%set(gca,'YTick',[])
ylab=ylabel({'Forward';'Particles'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlabel('Time (ms)')

%%
if numel(Sim.pfs)==2
    subplot(Nrows,Ncols,2), cla, hold on
    h=fill([Sim.tvec Sim.tvec(ind)],[R.p zeros(1,Sim.K)],gray);
    set(h,'edgecolor',gray)
    axis([min(Sim.tvec(xs)) max(Sim.tvec(xs)) 0 1])
    set(gca,'XTickLabel',[]), set(gca,'YTickLabel',[])
    set(gca,'XTick',[]), set(gca,'YTick',[])
    title('Prior Sampler')

    subplot(Nrows,2,4), cla, hold on
    plot(Sim.tvec,(S2.C'-cmin)/(cmax-cmin)+1,'Color',ccol(3,:))
    plot(Sim.tvec,(Cbar2-cmin)/(cmax-cmin)+1,'Color',col(3,:),'linewidth',2)
    plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',1,'markersize',4)
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')
    stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',2)
    BarVar2=Ibar2+Ivar2;
    BarVar2(BarVar2>1)=1;
    stem(Sim.tvec,BarVar2,'Marker','none','Color',ccol(3,:),'LineWidth',2)
    stem(Sim.tvec,Ibar2,'Marker','none','Color',col(3,:),'LineWidth',2)

    axis(AX)
    %set(gca,'XTickLabel',[]),
    set(gca,'YTickLabel',[])
    %set(gca,'XTick',[]),
    %set(gca,'YTick',[])
    xlabel('Time (ms)')
end

fig=figure(2);
bgr=0.5*[7 7];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior;