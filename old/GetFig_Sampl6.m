function GetFig_Sampl6(Sim,R,S,M)

%% preset stuff for fig
figure(2), cla, clf,
set(gcf, 'color', 'w');

gray    = [0.75 0.75 0.75];         %define gray
col     = [1 0 0; 0 .5 0];          %define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   %define colors for std

O           = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
ONaNind     = find(~isfinite(O));
Oind        = find(isfinite(O));
O(ONaNind)  = [];

nind    = find(R.n==1);
xs      = [Sim.tvec(nind(3)-20) Sim.tvec(nind(3)+20)];    %set the limits of the x-axis
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
l1   = .25;     %left of 1st col
w   = .35;      %width of subfigs
b   = .21;      %bottom of subfigs
h   = .3;       %height of subfigs
l2  = .05+l1+w; %left of second col
yfs = 9;       %ylabel font size
xfs = 9;       %xlabel font size
tfs = 12;       %title font size

%% make plots
for i=1:2
    if i==1,
        subplot(2,2,1), %subplot('Position',[l1 b w h]), %
        cla, hold on
        title('Forward-Only','fontsize',tfs),
        ylab=ylabel({'Prior'; 'Sampler'});
        set(ylab,...
             'fontsize',yfs, ...
             'Rotation',0,...
             'HorizontalAlignment','right',...
             'verticalalignment','middle')
    end
    if i==2,
        subplot(2,2,3),%subplot('Position',[l2 b w h]),
        cla, hold on
        ylab=ylabel({'Backwards';'Sampler'});
        set(ylab,...
             'fontsize',yfs, ...
             'Rotation',0,...
             'HorizontalAlignment','right',...
             'verticalalignment','middle')
        xlabel('Time (sec)', 'fontsize',xfs);
    end

    plot(Sim.tvec,(S{i}.C'-cmin)/(cmax-cmin)+1,'Color',ccol(i,:))                 %plot calcium particles
    plot(Sim.tvec,(fCbar(i,:)-cmin)/(cmax-cmin)+1,'Color',col(i,:),'linewidth',2) %plot calcium forward mean
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)              %plot true calcium

    plot(Oind*Sim.dt,(O-cmin)/(cmax-cmin)+1,'.k','LineWidth',1,'markersize',7)

    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    BarVar=fnbar(i,:)+fnvar(i,:);                                                   %make var of spikes not exceed 1
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',1)                   %plot true spikes
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(i,:),'LineWidth',2)         %plot forward spike var
    stem(Sim.tvec,fnbar(i,:),'Marker','none','Color',col(i,:),'LineWidth',2)      %plot forward spike mean
    axis(AX)
    set(gca,'YTickLabel',[])%,'XTickLabel',[])

    if i==1,
        subplot(2,2,2),
        title('Forward-Backward','fontsize',tfs),
    else
        subplot(2,2,4),
        xlabel('Time (sec)', 'fontsize',xfs);
    end,
    cla, hold on
    hfill=fill([Sim.tvec Sim.tvec(ind)],([M(i).Cbar-sqrt(M(i).Cvar) M(i).Cbar(ind)+sqrt(M(i).Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(i,:));
    set(hfill,'edgecolor',ccol(i,:))
    plot(Sim.tvec,(M(i).Cbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(i,:))
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
    plot(Oind*Sim.dt,(O-cmin)/(cmax-cmin)+1,'.k','LineWidth',1,'markersize',7)
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    BarVar=M(i).nbar+M(i).nvar;                                                   %make var of spikes not exceed 1
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',1)                   %plot true spikes
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(i,:),'LineWidth',2)         %plot forward spike var
    stem(Sim.tvec,M(i).nbar,'Marker','none','Color',col(i,:),'LineWidth',2)      %plot forward spike mean
    axis(AX)
    set(gca,'YTickLabel',[])%,'XTickLabel',[])
end

%% add some labels
% text('Interpreter','tex',...
%     'String',{'Calcium';'Paricles'},...
%     'Position',[-1.2 1.543 17.32],...
%     'FontSize',yfs,...
%     'FontName','Helvetica')
% 
% text('Interpreter','tex',...
%     'String',{'Spike';'Histogram'},...
%     'Position',[-1.2 0.4749 17.32],...
%     'FontSize',yfs,...
%     'FontName','Helvetica')

%%print to (color) eps
fig=figure(2);
wh=[7 4.5];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\sampl