function GetFig_Sampl8(Sim,R,S,M,n)

%% preset stuff for fig
figure(n), cla, clf,
%set(gcf, 'color', 'w');

gray    = [0.75 0.75 0.75];         %define gray
col     = [1 0 0; 0 .5 0];          %define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   %define colors for std

O           = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
ONaNind     = find(~isfinite(O));
Oind        = find(isfinite(O));
O(ONaNind)  = [];

nind    = find(R.n);             %find spike times
nind    = nind(n);
shift   = 20;
xs      = [Sim.tvec(nind(1)-shift) Sim.tvec(nind(1)+shift)];%set the limits of the x-axis
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
fs  = 15;       %default font size
yfs = fs;       %ylabel font size
xfs = fs;       %xlabel font size
titfs = fs+3;     %title font size
ticfs = fs;     %tick font size
texfs = fs;     %text font size
tl  = [.03 .25];%tick length
sw  = 5;        %spike width
xticks = AX(1):Sim.dt*Sim.freq*2:AX(2);
tx  = 1.435;
ty  = 2.1;

%% make plots
for i=1:2

    %% plot particles
    if i==1,                                                                        %top left plot
        subplot('Position',[l1 b2 w h]), %subplot(2,2,1)
        cla, hold on
        title('Particles','fontsize',titfs),
        ylab=ylabel({'Prior'; 'Sampler'});
        set(gca,'XTick',xticks,'XTickLabel',[],'TickLength',tl)
        %         text(tx,ty,'(A)','fontsize',texfs);
    end
    if i==2,                                                                        %top right plot
        subplot('Position',[l1 b w h]), %subplot(2,2,3)
        cla, hold on
        ylab=ylabel({'Conditional';'Sampler'});
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',xticks-AX(1),'TickLength',tl,'fontsize',ticfs)
        %         text(tx,ty,'(C)','fontsize',texfs);
    end
    set(ylab,'fontsize',yfs)%'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis(AX)
    set(gca,'YTickLabel',[])

    % calcium plot
    plot(Sim.tvec,(S{i}.C'-cmin)/(cmax-cmin)+1,'Color',ccol(i,:))                   %plot calcium particles
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',2)              %plot true calcium
    plot(Oind*Sim.dt,(O-cmin)/(cmax-cmin)+1,'.k','LineWidth',1,'markersize',7)      %plot observations
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    % spike plot
    stem(Sim.tvec,S{i}.n','Marker','none','Color',ccol(i,:),'LineWidth',1)                  %plot true spikes

    %% plot inferred distributions
    if i==1,                                                                        %bottom left
        subplot('Position',[l2 b2 w h]), %subplot(2,2,2),
        cla, hold on
        title('Inferred Distributions','fontsize',titfs),
        set(gca,'XTickLabel',[])
        set(gca,'XTick',xticks,'XTickLabel',[],'TickLength',tl)
        %         text(tx,ty,'(B)','fontsize',texfs);
    else                                                                            %bottom right
        subplot('Position',[l2 b w h]),%subplot(2,2,4),
        cla, hold on
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',xticks-AX(1),'TickLength',tl,'fontsize',ticfs)
        %         text(tx,ty,'(D)','fontsize',texfs);
    end
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
% fig=figure(2);
% wh=[7 4.5];
% set(fig,'PaperPosition',[0 11-wh(2) wh]);
% print -depsc C:\D\Research\liam\SMC_EM_GLM\sampl