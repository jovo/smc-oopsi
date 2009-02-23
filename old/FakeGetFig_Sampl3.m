function FakeGetFig_Sampl3(Sim,R,S,M)

%% preset stuff for fig
figure(2), cla, clf,
set(gcf, 'color', 'w');

gray    = [0.75 0.75 0.75];         %define gray
col     = [1 0 0; 0 .5 0];          %define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   %define colors for std

O       = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%make O only be observations (NaN'ing R.O at not observation times)
Oind    = find(~isnan(O));          %find the indices of observations that are not NaN's

xs      = [Sim.tvec(Oind(2)) Sim.tvec(Oind(end-2))];    %set the limits of the x-axis
xmin    = find(Sim.tvec>xs(1),1);   %index of lower limit on x
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;%index of upper limit on x
xind    = xmin:xmax;                %indices of x-axis
AX      = [xs 0 2];                 %axes

%get min and max of calcium to normalize within plots
for i=1:2
    cmin(i) = min(min(S(i).C(:,xind)));
    cmax(i) = max(max(S(i).C(:,xind)));
end
cmin=min(cmin(:));
cmax=max(cmax(:));

%get forward means and variances
for i=1:2
    fCbar(i,:) = sum(S(i).w_b.*S(i).C,1);
    fnbar(i,:) = sum(S(i).w_b.*S(i).n,1);
    fnvar(i,:) = sum((repmat(M(i).nbar,Sim.N,1)-S(i).n).^2)/Sim.N;
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
        subplot('Position',[l1 b w h]), cla, hold on
        title('Prior','fontsize',tfs), 
    end
    if i==2, 
        subplot('Position',[l2 b w h]), cla, hold on
        title('Backwards','fontsize',tfs),     
    end

    plot(Sim.tvec,(S(i).C'-cmin)/(cmax-cmin)+1,'Color',ccol(i,:))                 %plot calcium particles
    plot(Sim.tvec,(fCbar(i,:)-cmin)/(cmax-cmin)+1,'Color',col(i,:),'linewidth',2) %plot calcium forward mean
    plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',.5,'markersize',3)        %plot observations
    plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)              %plot true calcium
    plot(Sim.tvec,ones(size(Sim.tvec)),'k')

    BarVar=fnbar(i,:)+fnvar(i,:);                                                   %make var of spikes not exceed 1
    BarVar(BarVar>1)=1;
    stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',1)                   %plot true spikes
    stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(i,:),'LineWidth',2)         %plot forward spike var
    stem(Sim.tvec,fnbar(i,:),'Marker','none','Color',col(i,:),'LineWidth',2)      %plot forward spike mean
    axis(AX)
    set(gca,'YTickLabel',[])%,'XTickLabel',[])
    
    xlab=xlabel('Time (sec)');
    set(xlab,'fontsize',xfs)
end

%% add some labels
text('Interpreter','tex',...
	'String',{'Calcium';'Paricles'},...
	'Position',[-1.2 1.543 17.32],...
	'FontSize',yfs,...
    'FontName','Helvetica')
    
text('Interpreter','tex',...
	'String',{'Spike';'Histogram'},...
	'Position',[-1.2 0.4749 17.32],...
	'FontSize',yfs,...
    'FontName','Helvetica')

%%print to (color) eps
fig=figure(2);
wh=[6 2]
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\sampl