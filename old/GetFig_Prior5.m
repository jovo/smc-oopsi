function GetFig_Prior5(Sim,R,M,spt,priors)

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

b1  = .13;                      %bottom
h   = .225;                     %height
hs  = .05;                      %height space
b2  = b1+h+hs;                  %bottom of 2nd row
b3  = b2+h+hs;                  %bottom of 3rd row

%other stuff
tl  = [.04 .25];                %tick length [2d, 3d]
yfs = 14;                       %ylabel font size
xfs = 14;                       %xlabel font size
tfs = 19;                       %title font size
sw  = 2;                        %spike width

%% make plots
for i=1:3
    
    % plot prior
    subplot('Position',[l1+(w+ws)*(i-1) b3 w h])            %subplot('Position',[left bottom width height])  subplot(3,3,i) %
    Sim.x       = ones(Sim.T,1);                            %extize stimulus
    Sim.x(spt)  = priors(i);                                %add appropriate pulese
    stem(Sim.tvec,Sim.x,'color','k','marker','none','linewidth',sw)
    axis(priAX)
    set(gca,'YTickLabel',[],'XTickLabel',[],'box','off')
    if i==1                                                 %top left plot
        ylab=ylabel({'External';'Stimulus'});               %gets ylabel
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs,'color','k')
    end                                                     %the other prior figs get no ylabel
        set(gca,'YTick',priors,'YTickLabel',[],'TickLength',tl,'XTick',xticks)
    if i==1                                                 %label the columns
        title('Flat Prior','fontsize',tfs)
        set(gca, 'YTickLabel',priors)
    elseif i==2
        title('Weak Prior','fontsize',tfs)
    elseif i==3
        title('Strong Prior','fontsize',tfs)
    end

    % plot observations
    subplot('Position',[l1+(w+ws)*(i-1) b2 w h])            %subplot('Position',[left bottom width height]) subplot(3,3,i+Ncols) %
    cla, hold on
    plot(Oind,(O-omin)/(omax-omin),'.k','LineWidth',1,'markersize',10)
    set(gca,'YTick',[0:.5:1],'YTickLabel',[],'TickLength',tl,'box','off')
    axis(OAX)
    if i==1                                                 %left most plot gets ylabel
        ylab=ylabel('Observations');
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs,'color','k')
    end
    set(gca,'XTick',[0:40:240],'XTickLabel',[])

    % plot spike inference
    subplot('Position',[l1+(w+ws)*(i-1) b1 w h])            %subplot(3,3,i+Ncols*2)
    cla, hold on
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
    BarVar=M(i).nbar+M(i).nvar;
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',sw)
    stem(Sim.tvec,M(i).nbar,'Marker','none','Color',col(2,:),'LineWidth',sw)
    axis([parAX])
    
    set(gca,'YTick',[0 .5 1], 'YTickLabel',[],'XTick',xticks,'TickLength',tl)
    if i==1
        ylab=ylabel({'Inferred';'Spikes'});
        set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(2,:),'fontsize',yfs)
    end
    xlab=xlabel('Time (sec)');
    set(xlab,'fontsize',xfs);
end

%% print to (color) eps
fig=figure(4);
wh=[7 3];   %width and height
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\prior