function GetFig_Sampl7(Sim,R,S,M)

%% preset stuff for fig
figure(2), cla, clf,
%set(gcf, 'color', 'w');

gray    = [0.75 0.75 0.75];         %define gray
col     = [1 0 0; 0 .5 0];          %define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   %define colors for std

O           = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
ONaNind     = find(~isfinite(O));
Oind        = find(isfinite(O));
O(ONaNind)  = [];

nind    = find(R.n==1);             %find spike times
nind    = nind(5);
xs      = [Sim.tvec(nind-20) Sim.tvec(nind+20)];%set the limits of the x-axis
xmin    = find(Sim.tvec>xs(1),1);   %index of lower limit on x
xmax    = find(Sim.tvec>xs(2),1)-1; %index of upper limit on x
xind    = xmin:xmax;                %indices of x-axis
AX      = [xs 0 2];                 %axes
ind     = Sim.T:-1:1;

%get min and max of calcium to normalize within plots
for i=1:2
    cmin(i) = min(min(M(i).Cbar(xind)-sqrt(M(i).Cvar(xind))));
    cmax(i) = max(max(M(i).Cbar(xind)+sqrt(M(i).Cvar(xind))));
end
cmin=min(cmin(:));
cmax=max(cmax(:));

%get forward means and variances
for i=1:2
    fCbar(i,:) = sum(S{i}.w_f.*S{i}.C,1);
    fnbar(i,:) = sum(S{i}.w_f.*S{i}.n,1);
    fnvar(i,:) = sum((repmat(fnbar(i,:),Sim.N,1)-S{i}.n).^2)/Sim.N;
end

% set subfig sizes
l1  = .1;      %left of 1st col
w   = .41;      %width of subfigs
b   = .22;      %bottom of subfigs
h   = .3;       %height of subfigs
b2  = b+h+.05;  %bottom of 2nd row
l2  = .05+l1+w; %left of second col
yfs = 13;       %ylabel font size
xfs = 13;       %xlabel font size
tfs = 13;       %title font size
tl  = [.05 .25];%tick length
sw  = 5;        %spike width
xticks = AX(1):Sim.dt*Sim.freq*2:AX(2);
tx  = 1.435;
ty  = 2.1;

%% make plots
for i=1:2
    
    %% plot forward only stuff
    if i==1,                                                                        %top left plot
        subplot('Position',[l1 b2 w h]), %subplot(2,2,1)
        cla, hold on
        title('Forward-Only','fontsize',tfs),
        ylab=ylabel({'Prior'; 'Sampler'});
        set(gca,'XTick',xticks,'XTickLabel',[],'TickLength',tl)
        text(tx,ty,'(A)','fontsize',tfs);
    end
    if i==2,                                                                        %top right plot
        subplot('Position',[l1 b w h]), %subplot(2,2,3)
        cla, hold on
        ylab=ylabel({'Conditional';'Sampler'});
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',xticks-AX(1),'TickLength',tl)
        text(tx,ty,'(C)','fontsize',tfs);
    end
    set(ylab,'fontsize',yfs)%'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis(AX)
    set(gca,'YTickLabel',[])

    % calcium plot
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',2)              %plot true calcium
    plot(Sim.tvec,(S{i}.C'-cmin)/(cmax-cmin)+1,'Color',ccol(i,:))                   %plot calcium particles
    plot(Sim.tvec,(fCbar(i,:)-cmin)/(cmax-cmin)+1,'Color',col(i,:),'linewidth',2)   %plot calcium forward mean
    plot(Oind*Sim.dt,(O-cmin)/(cmax-cmin)+1,'.k','LineWidth',1,'markersize',7)      %plot observations
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    % spike plot
    BarVar=fnbar(i,:)+fnvar(i,:);                                                   %make var of spikes not exceed 1
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)                  %plot true spikes
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(i,:),'LineWidth',sw)          %plot forward spike var
    stem(Sim.tvec,fnbar(i,:),'Marker','none','Color',col(i,:),'LineWidth',sw)       %plot forward spike mean

    %% plot forward-backward stuff
    if i==1,                                                                        %bottom left
        subplot('Position',[l2 b2 w h]), %subplot(2,2,2),
        cla, hold on
        title('Forward-Backward','fontsize',tfs),
        set(gca,'XTickLabel',[])
        set(gca,'XTick',xticks,'XTickLabel',[],'TickLength',tl)
        text(tx,ty,'(B)','fontsize',tfs);
    else                                                                            %bottom right
        subplot('Position',[l2 b w h]),%subplot(2,2,4),
        cla, hold on
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',xticks-AX(1),'TickLength',tl)
        text(tx,ty,'(D)','fontsize',tfs);
    end,
    axis(AX)
    set(gca,'YTickLabel',[])

    % calcium plot
    hfill=fill([Sim.tvec Sim.tvec(ind)],([M(i).Cbar-sqrt(M(i).Cvar) M(i).Cbar(ind)+sqrt(M(i).Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(i,:));
    set(hfill,'edgecolor',ccol(i,:))
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',2)
    plot(Sim.tvec,(M(i).Cbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(i,:))
    plot(Oind*Sim.dt,(O-cmin)/(cmax-cmin)+1,'.k','LineWidth',1,'markersize',7)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    % spike plot
    BarVar=M(i).nbar+M(i).nvar;                                                     %make var of spikes not exceed 1
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)                  %plot true spikes
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(i,:),'LineWidth',sw)          %plot forward spike var
    stem(Sim.tvec,M(i).nbar,'Marker','none','Color',col(i,:),'LineWidth',sw)        %plot forward spike mean
end


%% print to (color) eps
fig=figure(2);
wh=[7 4.5];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\bernoulli\sampl