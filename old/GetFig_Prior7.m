function GetFig_Prior7(Sim,R,M,spt,priors)

%% preset figure stuff
figure(4), clf,
set(gcf, 'color', 'w');

%let O be only observations at sample times
O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
Onan    = find(~isfinite(O));
Oind    = find(isfinite(O));
O(Onan) = [];

%limits of x-axis
xs      = [0 Sim.Nsec-2*(Sim.dt*Sim.freq)];
xticks  = Sim.tvec(1)+Sim.dt:Sim.dt*40:Sim.tvec(end);
xticks  = xticks-xticks(1);

%axes
priAX   = [xs 0 max(priors)];   %axes for prior plots
parAX   = [xs 0 1];             %axes for particle plots
OAX     = [xs/Sim.dt 0 1];      %axes for observation plots

%get min and max of calcium for normalizing the plots
omin    = min(O);
omax    = max(O);

%define colors
gray=[0.75 0.75 0.75];          %define gray
col=[0 0 1; 0 .5 0];            %define colors for mean
ccol=col+.8; ccol(ccol>1)=1;    %define colors for std

%subfig positions
l1   = .225;                    %left side of fig
w   = .21;                      %width
ws  = .05;                      %width spacer

b   = .15;                      %bottom
h   = .15;                     %height
hs  = .05;                      %height space

%other stuff
tl  = [.04 .25];                %tick length [2d, 3d]
yfs = 15;                       %ylabel font size
xfs = 15;                       %xlabel font size
tfs = 15;                       %title font size
sw  = 2;                        %spike width

%% make plots
for i=1:3

    % plot prior
    subplot('Position',[l1+(w+ws)*(i-1) b+(h+hs)*3 w h])            %subplot('Position',[left bottom width height])  subplot(3,3,i) %
    Sim.x       = ones(Sim.T,1);                            %extize stimulus
    if i==2, Sim.x(spt(1):spt(end)) = 5*sin(pi:pi/(spt(end)-spt(1)):2*pi)+5;                                %add appropriate pulese
    elseif i==3, Sim.x(spt)         = priors(i);
    end
    stem(Sim.tvec,Sim.x,'color','k','marker','none','linewidth',sw)
    axis(priAX)
    if i==1                                                 %top left plot
        ylab=ylabel({'External';'Stimulus'});               %gets ylabel
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs,'color','k')
    end                                                     %the other prior figs get no ylabel
    set(gca,'YTick',[1 3 5],'YTickLabel',[],...
            'XTickLabel',[],'XTick',xticks,...
            'box','off','TickLength',tl)
    if i==1                                                 %label the columns
        title({'None'},'fontsize',tfs)
    elseif i==2
        title({'Coarse'},'fontsize',tfs)
    elseif i==3
        title({'Fine'},'fontsize',tfs)
    end

    % plot observations
    subplot('Position',[l1+(w+ws)*(i-1) b+(h+hs)*2 w h])            %subplot('Position',[left bottom width height]) subplot(3,3,i+Ncols) %
    cla, hold on
    plot(Oind,(O-omin)/(omax-omin),'.k','LineWidth',1,'markersize',10)
    set(gca,'YTick',[0:.5:1],'YTickLabel',[],'TickLength',tl,'box','off')
    axis(OAX)
    if i==1                                                 %left most plot gets ylabel
        ylab=ylabel('Observations');
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs,'color','k')
    end
    set(gca,'XTick',[0:40:240],'XTickLabel',[])

    % plot with spike history terms
    m=2;
    subplot('Position',[l1+(w+ws)*(i-1) b+(h+hs) w h])            %subplot(3,3,i+Ncols*2)
    cla, hold on
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
    BarVar=M(m,i).nbar+M(m,i).nvar;
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',sw)
    stem(Sim.tvec,M(m,i).nbar,'Marker','none','Color',col(2,:),'LineWidth',sw)
    axis([parAX])
    set(gca,'YTick',[0 .5 1], 'YTickLabel',[],'XTick',xticks,'XTickLabel',[],'TickLength',tl)
    if i==1
        ylab=ylabel({'With';'Spike Histories'});
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(2,:),'fontsize',yfs)
    end
    
    % plot without spike history terms
    m=1;
    subplot('Position',[l1+(w+ws)*(i-1) b w h])            %subplot(3,3,i+Ncols*2)
    cla, hold on
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
    BarVar=M(m,i).nbar+M(m,i).nvar;
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(1,:),'LineWidth',sw)
    stem(Sim.tvec,M(m,i).nbar,'Marker','none','Color',col(1,:),'LineWidth',sw)
    axis([parAX])
    set(gca,'YTick',[0 .5 1], 'YTickLabel',[],'XTick',xticks,'TickLength',tl)
    if i==1
        ylab=ylabel({'Without';'Spike Histories'});
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(1,:),'fontsize',yfs)
    end
    
    xlab=xlabel('Time (sec)');
    set(xlab,'fontsize',xfs);
end

%text('String','Knowledge',...
%    'Position',[-2.1 5.4 17.32],...
%    'FontSize',tfs);


%% print to (color) eps
fig=figure(4);
wh=[7 3];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior