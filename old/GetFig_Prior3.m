function GetFig_Prior3(Sim,R,S,M,spt,priors)

%% preset figure stuff
figure(4), clf,
set(gcf, 'color', 'w');

O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
Onan    = O;
Oind    = find(~isfinite(O));
O(Oind) = [];
Ointer  = interp1(0:Sim.freq:Sim.T,[0 O],Sim.tvec/Sim.dt,'linear');

xs      = [0 Sim.Nsec-2*(Sim.dt*Sim.freq)];                                          %limits of x-axis
xmin    = find(Sim.tvec>xs(1),1);                           %min of x-axis
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;                    %max of x-axis
xind    = xmin:xmax;                                        %indices of limits of x-axis

priAX   = [xs 0 max(priors)];                                         %axes for prior plots
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
    Sim.x       = ones(Sim.T,1);                %extize stimulus
    Sim.x(spt)  = priors(i);         %add appropriate pulese
    stem(Sim.tvec,Sim.x,'color','k','marker','none','linewidth',4)
    axis(priAX)
    set(gca,'YTickLabel',[],'XTickLabel',[],'box','off')
    if i==1
        ylab=ylabel({'External';'Stimulus'});
        set(ylab,...
            'Rotation',0,...
            'HorizontalAlignment','right',...
            'verticalalignment','middle',...
            'fontsize',yfs,...
            'color','k')
    end
    set(gca,...
        'YTick',priors,...
        'YTickLabel',[]...
        ,'TickLength',tl...
        ,'XTick',[0:.1:.6])
    if i==1
        title('Flat Prior','fontsize',tfs)
        set(gca, 'YTickLabel',priors)
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
    plot(Sim.tvec,(Ointer-cmin)/(cmax-cmin)+1,'-k','LineWidth',1,'markersize',10)
    plot(Sim.tvec,(Onan-cmin)/(cmax-cmin)+1,'.k','LineWidth',1,'markersize',10)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',1)
    BarVar=M(i).nbar+M(i).nvar;
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',3)
    stem(Sim.tvec,M(i).nbar,'Marker','none','Color',col(2,:),'LineWidth',3)

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