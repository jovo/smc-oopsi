n=1;
figure(2), cla, clf,

gray    = [0.75 0.75 0.75];         %define gray
col     = [1 0 0; 0 .5 0];          %define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   %define colors for std

O           = R.F.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
ONaNind     = find(~isfinite(O));
Oind        = find(isfinite(O));
O(ONaNind)  = [];
finv        = ((P.k_d.*(P.beta-O))./(O-P.beta-P.alpha)).^(1/P.n);

nind    = find(R.n);                %find spike times
nind    = nind(n);
xmin    = Oind(find(Oind>nind,1)-1)-1;
xmax    = Oind(find(Oind>nind,1))+1;
xs      = Sim.tvec([xmin xmax]);%set the limits of the x-axis
xind    = xmin:xmax;                %indices of x-axis
ind     = Sim.T:-1:1;               %inverse index for 'fill' plots

%get min and max of calcium to normalize within plots
cshift=inf;
for i=1:2
    cmin(i) = min(min(R.C(xind)),min(min(S{i}.C(:,xind))));
    cmax(i) = max(max(R.C(xind)),max(max(S{i}.C(:,xind))));
    hmin(i) = min(min(P.omega*R.h(xind)),min(min(P.omega*S{i}.h(:,xind))));
    hmax(i) = max(max(P.omega*R.h(xind)),max(max(P.omega*S{i}.h(:,xind))));
    cshift  = min(min(min(S{i}.C(:,xind))),cshift);
end
cmin    = min(cmin(:));
cmax    = max(cmax(:));
cdiff   = cmax-cmin;

hmin    = min(hmin(:));
hmax    = max(hmax(:));
hdiff   = hmax-hmin;

for i=1:2
   C{i}=(S{i}.C(:,xind)-cmin)/cdiff;
   hh{i}=(P.omega*S{i}.h(:,xind)-hmin)/hdiff;
end
    

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
titfs = fs+3;   %title font size
ticfs = fs;     %tick font size
texfs = fs;     %text font size
tl  = [.03 .03];%tick length
sw  = 5;        %spike width
bw  = .3;
bw2 = .0015;
xticks = xs(1):Sim.dt:xs(2);
tx  = 1.435;
ty  = 2.1;
sp  = 1.2;
lw  = 2;

AX      = [xs 0 3*sp];                 %axes
% [foo ind]=find(R.n==0);
% if ~isempty(ind), R.n(ind)=NaN; end

% make plots
for i=1:2

    %% plot particles
    if i==1,                                                                        %top left plot
        subplot('Position',[l1 b2 w h]), %subplot(2,2,1)
        cla, hold on
        title('Particles','fontsize',titfs),
        ylab=ylabel({'Prior'; 'Sampler'});
        set(gca,'XTick',xticks,'XTickLabel',[],'TickLength',tl,'YTick',[])
    end
    if i==2,                                                                        %top right plot
        subplot('Position',[l1 b w h]), %subplot(2,2,3)
        cla, hold on
        ylab=ylabel({'Conditional'; 'Sampler'});
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',[{'u-1'}, {'u'}, {'u+1'}, {'u+2'}, {'v'}, {''}],'TickLength',tl,'fontsize',ticfs,'YTick',[])        %         text(tx,ty,'(C)','fontsize',texfs);
    end
    set(ylab,'fontsize',yfs)%'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTickLabel',[])
    
    % spike particles
    for x=xind(1)+1:xind(end)
    	bar(Sim.tvec(x),R.n(x)+2*sp,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw2)%plot true spikes
        bar(Sim.tvec(x),sum(S{i}.n(:,x))/5+2*sp,'EdgeColor',col(i,:),'FaceColor',col(i,:),'BarWidth',bw2)                  %plot spike particles
      	bar(Sim.tvec(x),2*sp,'EdgeColor','w','FaceColor','w','BarWidth',.003)%plot true spikes
    end
    plot(Sim.tvec,0*ones(size(Sim.tvec)),'k','linewidth',lw)
    plot(Sim.tvec,sp*ones(size(Sim.tvec)),'k','linewidth',lw)
    plot(Sim.tvec,2*sp*ones(size(Sim.tvec)),'k','linewidth',lw)

    % calcium particles
    plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2)              %plot true calcium
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),C{i}(nn,x)','.','Color',col(i,:),'markersize',10*exp(S{i}.w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),C{i}(:,1:length(xind))','Color',col(i,:))
%     plot(Oind*Sim.dt,finv-cshift,'.k','LineWidth',1,'markersize',7)      %plot observations

    % h particles
    plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff+sp,'color',gray,'LineWidth',2)              %plot true calcium
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),hh{i}(nn,x)'+sp,'.','Color',col(i,:),'markersize',10*exp(S{i}.w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),hh{i}(:,1:length(xind))'+sp,'Color',col(i,:))
    axis(AX)

    
    %% plot inferred distributions
    if i==1,                                                                        %bottom left
        subplot('Position',[l2 b2 w h]), %subplot(2,2,2),
        cla, hold on
        title('Inferred Distributions','fontsize',titfs),
        set(gca,'XTickLabel',[])
        set(gca,'XTick',xticks,'XTickLabel',[],'TickLength',tl,'YTick',[])
    else                                                                            %bottom right
        subplot('Position',[l2 b w h]),%subplot(2,2,4),
        cla, hold on
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',[{'u-1'}, {'u'}, {'u+1'}, {'u+2'}, {'v'}, {'v+1'}],'TickLength',tl,'fontsize',ticfs,'YTick',[])        %         text(tx,ty,'(C)','fontsize',texfs);
    end
    axis(AX)
    set(gca,'YTickLabel',[])

    % calcium plot
    ptiles = GetPercentiles([.05 .95],S{i}.w_b,S{i}.C);
    plot(Sim.tvec,R.C+cdiff-cshift,'color',gray,'LineWidth',2)
    hfill=fill([Sim.tvec Sim.tvec(ind)],[ptiles(1,:) ptiles(2,ind)]+cdiff-cshift,ccol(i,:));
    set(hfill,'edgecolor',ccol(i,:))
    plot(Sim.tvec,M(i).Cbar+cdiff-cshift,'linewidth',2,'color',col(i,:))
    plot(Oind*Sim.dt,finv+cdiff-cshift,'.k','LineWidth',1,'markersize',7)      %plot observations
    plot(Sim.tvec,cdiff*ones(size(Sim.tvec)),'k')

    % spike plot
    BarVar=M(i).nbar+M(i).nvar;                                                     %make var of spikes not exceed 1
    BarVar(BarVar>1)=1;
	bar(Sim.tvec,R.n*cdiff,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
	bar(Sim.tvec,BarVar*cdiff,'EdgeColor',ccol(i,:),'FaceColor',ccol(i,:),'BarWidth',bw)%plot true spikes
	bar(Sim.tvec,M(i).nbar*cdiff,'EdgeColor',col(i,:),'FaceColor',col(i,:),'BarWidth',bw)%plot true spikes
    plot(Sim.tvec,0*ones(size(Sim.tvec)),'k')
end

% print to (color) eps
fig=figure(2);
wh=[7 4.5];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\sampl3