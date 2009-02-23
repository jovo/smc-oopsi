function GetFig_Prior2(Sim,R,S,M)

%% preset figure stuff
figure(4), clf,
set(gcf, 'color', 'w');

O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
xs      = [0 .65];                                          %limits of x-axis
xmin    = find(Sim.tvec>xs(1),1);                           %min of x-axis
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;                    %max of x-axis
xind    = xmin:xmax;                                        %indices of limits of x-axis

priAX   = [xs 0 1];                                         %axes for prior plots
parAX   = [xs 0 2];                                         %axes for particle plots

%get min and max of calcium for normalizing the ploots
cmin   = min(min(M(1).Cbar(:,xind)-sqrt(M(1).Cvar(:,xind))));
cmax   = max(max(M(1).Cbar(:,xind)+sqrt(M(1).Cvar(:,xind))));

%define colors
gray=[0.75 0.75 0.75];                                      %define gray
col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];   %define colors for mean
ccol=col+.8; ccol(ccol>1)=1;                                %define colors for std
ind=Sim.T:-1:1;                                             %inverse indices for 'fill' function


Nrows = 2;
Ncols = 3;

l   = .225;     %left side of fig
b   = .15;      %bottom
w   = .21;      %width
h   = .28;      %height
h2  = .15;      %height of top row
ws  = .05;      %width spacer
hs  = .05;      %height space
tl  = [.04 .25];%tick length [2d, 3d]
yfs = 16;       %ylabel font size
xfs = 14;       %xlabel font size
tfs = 19;       %title font size

%% make plots
for i=1:3
    subplot('Position',[l+(w+ws)*(i-1) b+h+hs w h2]) %subplot('Position',[left bottom width height])
    stem(Sim.tvec,S(i).p(1,:)*4,'color',gray,'marker','none','linewidth',2)
    axis(priAX)
    set(gca,'YTickLabel',[],'XTickLabel',[],'box','off')
    if i==1
        ylab=ylabel({'Rate';'(Hz)'});
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs,'color',gray)
    end
    set(gca,'YTick',[.05 .15]*4,'YTickLabel',[],'TickLength',tl,'XTick',[0:.2:.6])    
    if i==1
        title('Flat Prior','fontsize',tfs)
        set(gca, 'YTickLabel',{'10';'30'})    
    elseif i==2
        title('Weak Prior','fontsize',tfs)
    elseif i==3
        title('Strong Prior','fontsize',tfs)
    end

    subplot('Position',[l+(w+ws)*(i-1) b w h]) %subplot('Position',[left bottom width height])
    cla, hold on
%     fillh=fill([Sim.tvec Sim.tvec(ind)],([M(i).Cbar-sqrt(M(i).Cvar) M(i).Cbar(ind)+sqrt(M(i).Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
%     set(fillh,'edgecolor',ccol(2,:))
%     plot(Sim.tvec,(M(i).Cbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))
%     plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
    plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'.k','LineWidth',1,'markersize',10)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',1)
    BarVar=M(i).nbar+M(i).nvar;
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',2)
    stem(Sim.tvec,M(i).nbar,'Marker','none','Color',col(2,:),'LineWidth',2)

    set(gca,'YTick',[], 'YTickLabel',[],'XTick',[0:.2:.6],'TickLength',tl)
    axis(parAX)
    if i==1
        ylab=ylabel({'Calcium';'and';'Spike Train'});
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(2,:),'fontsize',yfs)
    end
    xlab=xlabel('Time (sec)');
    set(xlab,'fontsize',xfs);
end

%% print to (color) eps
fig=figure(4);
wh=[6 3];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior