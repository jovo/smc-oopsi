function GetPriorFig4A(Sim,R,S,M)

O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.K_o);
Oind    = find(~isnan(O));
xs      = [0 .65];%[Sim.tvec(Oind(2)) Sim.tvec(Oind(end-1))];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];
priAX   = [xs 0 1];
parAX   = [xs 0 2];

cmin   = min(min(M(1).bCbar(:,xind)-sqrt(M(1).bCvar(:,xind))));
cmax   = max(max(M(1).bCbar(:,xind)+sqrt(M(1).bCvar(:,xind))));

gray    = [0.75 0.75 0.75];

col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8;
ccol(ccol>1)=1;
ind     = Sim.K:-1:1;

Nrows   = 2;
Ncols   = 3;

figure(4), clf,
set(gcf, 'color', 'w');

l1   = .25;
w   = .35;
b   = .21;
h   = .3;
hs  = .05;
l2  = .05+l1+w;
tl  = .05;
yfs = 16;
xfs = 14;
tfs = 19;

for i=1:3
    % subplot('Position',[l1 (b+h+hs) w h]) %subplot('Position',[left bottom width height])
    subplot(Nrows,Ncols,i)
    stem(Sim.tvec,S(i).p(1,:),'color',gray,'marker','none','linewidth',2)
    axis(priAX)
    set(gca,'YTickLabel',[])
    set(gca,'XTickLabel',[])
    if i==1
        ylab=ylabel({'Probability';'of Spiking'});
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs)
    end
    if i==1
        title('Flat Prior','fontsize',tfs)
    elseif i==2
        title('Weak Prior','fontsize',tfs)
    elseif i==3
        title('Strong Prior','fontsize',tfs)
    end

    % subplot('Position',[l1 b w h]) %subplot('Position',[left bottom width height])
    subplot(Nrows,Ncols,i+Ncols)
    cla, hold on
    h=fill([Sim.tvec Sim.tvec(ind)],([M(i).bCbar-sqrt(M(i).bCvar) M(i).bCbar(ind)+sqrt(M(i).bCvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
    set(h,'edgecolor',ccol(2,:))
    plot(Sim.tvec,(M(i).bCbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
    plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',1,'markersize',4)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
    BarVar=M(i).bIbar+M(i).bIvar;
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',2)
    stem(Sim.tvec,M(i).bIbar,'Marker','none','Color',col(2,:),'LineWidth',2)

    set(gca,'YTick',[]), set(gca,'YTickLabel',[])
    %     set(gca,'XTick',Sim.tvec(Oind)),
    axis(parAX)
    if i==1
        ylab=ylabel({'Calcium';'and';'Spike Train'});
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(2,:),'fontsize',yfs)
    end
    xlab=xlabel('Time (sec)');
    set(xlab,'fontsize',xfs);
end

fig=figure(4);
wh=[6 3];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior