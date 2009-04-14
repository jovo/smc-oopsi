% this m-file generates the data and then plots the fig demonstrating how
% the two different sampling strategies differ.  It generates the following:
%
% Sim:  simulation parameters
% P:    parameters of "real" neuron
% R:    "real" neuron data                  (smc_em_bern_real_exp)
% S:    simulation states for both samplers (smc_em_bern_main)
% M:    moments for both samplers           (smc_em_bern_main)
% fig:  see fig file for details            (GetSamplFig1)

%% start function
clear; clc;

[Sim P] = InitializeStuff;

% set simulation parameters
Sim.freq    = 3;                        %frequency of observations
Sim.Nsec    = 2.5;                      %# of sec
Sim.StimDim = 1;                        %# of stimulus dimensions
Sim.M       = 1;                        %number of spike history terms
Sim.N       = 5;
Sim.pf      = 1;                        %not vanilla particle filtering
Sim.T       = round(Sim.Nsec/Sim.dt);   %total # of steps (round deals with numerical error)
rem         = mod(Sim.T,Sim.freq);      %remainder
if rem~=0
    Sim.T=Sim.T-rem;                    %fix number of steps
end
Sim.T_o     = round(Sim.T/Sim.freq);    %number of observations (round deals with numerical error)
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec-Sim.dt*rem;%time vector
Sim.x       = ones(1,Sim.T);     %make one input stationary

P.k         = 0;                        %bias term
P.gamma     = .1e-5;                     %var gainh
P.zeta      = .1e-5;                    %var offset
P.A         = .08;
P.sigma_c   = .1;
P.tau_h     = .01;                %decay rate for spike history terms
P.a         = Sim.dt/P.tau_c;

% get "real" data
R.n         = zeros(1,Sim.T);   %spike times
R.h         = zeros(1,Sim.T);   %spike times
R.C         = P.C_init*ones(1,Sim.T);   %initialize calcium
epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
% spt         = [3*Sim.freq:Sim.freq:50 81 131 181 231 281 331 381 431 481];    %forced spike times
spt         = [181 231 281 331 381 431 481];    %forced spike times
spt         = spt(spt<Sim.T);
R.n(spt)    = 1;                %force spikes
epsilon_h = repmat(P.sigma_h*sqrt(Sim.dt),1,Sim.T).*randn(Sim.M,Sim.T); %generate noise on spike history

for t=2:Sim.T                   %update calcium
    R.h(:,t)= (1-Sim.dt./P.tau_h).*R.h(:,t-1)+R.n(t-1) + epsilon_h(:,t);%update h terms
    R.C(t)  = (1-P.a)*R.C(t-1) + P.A*R.n(t) + P.a*P.C_0 + epsilon_c(t);
end
F_mu        = P.alpha*Hill_v1(P,R.C)+P.beta;        %compute E[F_t]
F_var       = P.gamma*Hill_v1(P,R.C)+P.zeta;    %compute V[F_t]
R.F         = F_mu+sqrt(F_var).*randn(1,Sim.T);%add noise to observations
R.F(R.F<0)  = eps;                      %observations must be non-negative

% do EM recursion
for i=1:2
    if i==1, Sim.pf=0; else Sim.pf=1; end
    [S{i} M(i)] = smc_em_bern_FoBaMo_v5(Sim,R,P);
end

% figure(4), clf, hold on,
% plot(S{1}.C','b'),
% plot(R.C,'k','linewidth',2),
% plot(M(1).Cbar,'r','linewidth',2)

%% make a fig

n=1;
figure(2), cla, clf,

gray    = [0.75 0.75 0.75];         %define gray
col     = [1 0 0; .2 .2 1];          %define colors for mean
ccol    = col+.4; ccol(ccol>1)=1;   %define colors for std

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
    M(i).hbar = sum(S{i}.w_b.*S{i}.h,1);
    M(i).hvar = sum((repmat(M(i).hbar,Sim.N,1)-S{i}.h).^2)/Sim.N;
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
% w   = .41;      %width of subfigs
% b   = .22;      %bottom of subfigs
% h   = .3;       %height of subfigs
% b2  = b+h+.05;  %bottom of 2nd row
% l2  = .05+l1+w; %left of second col
fs  = 12;       %default font size
yfs = fs;       %ylabel font size
xfs = fs;       %xlabel font size
titfs = 12;   %title font size
ticfs = fs;     %tick font size
texfs = fs;     %text font size
tfs   = 14;
tl  = [.03 .03];%tick length
sw  = 5;        %spike width
bw  = .3;
bw2 = .0015;
xticks = xs(1):Sim.dt:xs(2);
tx  = 1.435;
ty  = 2.1;
sp  = 1.2;
lw  = 2;

Nrows=7;
Ncols=2;
AX  = [xs 0 1];
interp = 'none';

r2 = .5;
h  = .1;
l1 = .25;      %left of 1st col
hs = .05;
w  = .35;
l2 = l1+w+hs;
ms = 7;         %markersize
ms2= 18;
yshift = .05;

Sim.n = R.n;
Sim.n(Sim.n==0) = NaN;
obs = xind([2 5]);
for i=1:2

%     S{i}.n(S{i}.n<1e-3) = 0;
    % h particles
    subplot('Position',[l1 .34+(2-i)*r2 w h]) %[left bottom width height]    %     subplot(Nrows,Ncols,ii(1)),
    cla, hold on
    plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff,'color',gray,'LineWidth',2,'LineStyle','--')              %plot true calcium
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),hh{i}(nn,x)','.','Color',col(i,:),'markersize',50*(S{i}.w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),hh{i}(:,1:length(xind))','Color',col(i,:))
    set(gca,'YTick',[0 1],'YTickLabel',[]); %,round([hmin hmax]*100)/100)
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    title('Particles','fontsize',titfs),
    ylab=ylabel([{'Weighted'}; {'Spike History'}],'fontsize',yfs,'Interpreter',interp);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

    % spike particles
    subplot('Position',[l1 .22+(2-i)*r2 w h]) %[left bottom width height]     subplot(Nrows,Ncols,ii(2)),
    cla, hold on
    for x=xind(1)+1:xind(end)
        stem(Sim.tvec(x),Sim.n(x)-yshift,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
        bar(Sim.tvec(x),sum(S{i}.n(:,x))/5,'EdgeColor',col(i,:),'FaceColor',col(i,:),'BarWidth',bw2)                  %plot spike particles
    end
    set(gca,'YTick',[0 1],'YTickLabel',[]); %[0 1])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    ylab=ylabel([{'Spike Train'}],'fontsize',yfs,'Interpreter',interp);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

    % calcium particles
    subplot('Position',[l1 .1+(2-i)*r2 w h]) %[left bottom width height]    % subplot(Nrows,Ncols,ii(3)), cla, hold on
    hold on
    plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2,'LineStyle','--')              %plot true calcium
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),C{i}(nn,x)','.','Color',col(i,:),'markersize',50*(S{i}.w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),C{i}(:,1:length(xind))','Color',col(i,:))
%     plot(Sim.tvec(obs),(R.C(obs)-cmin)/cdiff,'.k','Markersize',ms2)
    set(gca,'YTick',[0 1],'YTickLabel',[]); %,round((([cmin cmax])-cmin)*100)/100)
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
    xlabel('Time (sec)', 'fontsize',xfs);
    ylab=ylabel([{'Calcium'}; {'Concentration'}],'fontsize',yfs,'Interpreter',interp);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis(AX)
    box off


    % h dist
    subplot('Position',[l2 .34+(2-i)*r2 w h]) %[left bottom width height]%    subplot(Nrows,Ncols,ii(1)+1),
    cla, hold on
    ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.h);
    hfill=fill([Sim.tvec Sim.tvec(ind)],(P.omega*[ptiles(1,:) ptiles(2,ind)]-hmin)/hdiff,ccol(i,:));
    set(hfill,'edgecolor',ccol(i,:))
    plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff,'color',gray,'LineWidth',2,'LineStyle','--')              %plot true calcium
    plot(Sim.tvec,(P.omega*M(i).hbar-hmin)/hdiff,'linewidth',2,'color',col(i,:))
    set(gca,'YTick',[0 1],'YTickLabel',[])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    title('Inferred Distributions','fontsize',titfs),


    % spike plot
    subplot('Position',[l2 .22+(2-i)*r2 w h]) %[left bottom width height]    %subplot(Nrows,Ncols,ii(2)+1),
    cla, hold on
    ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.n);
    %     bar(Sim.tvec,R.n,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
    bar(Sim.tvec,ptiles(2,:),'EdgeColor',ccol(i,:),'FaceColor',ccol(i,:),'BarWidth',bw)%plot true spikes
    bar(Sim.tvec,M(i).nbar,'EdgeColor',col(i,:),'FaceColor',col(i,:),'BarWidth',bw)%plot true spikes
    stem(Sim.tvec,Sim.n-yshift,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
    set(gca,'YTick',[0 1],'YTickLabel',[])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)


    % calcium plot
    subplot('Position',[l2 .1+(2-i)*r2 w h]) %[left bottom width height]    %subplot(Nrows,Ncols,ii(3)+1),
    cla, hold on
    ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.C);
    hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-cmin)/cdiff,ccol(i,:));
    set(hfill,'edgecolor',ccol(i,:))
    plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2,'LineStyle','--')
    plot(Sim.tvec,(M(i).Cbar-cmin)/cdiff,'linewidth',2,'color',col(i,:))
    plot(Sim.tvec(obs),(R.C(obs)-cmin)/cdiff,'.k','Markersize',ms2)
    set(gca,'YTick',[0 1],'YTickLabel',[])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
    xlabel('Time (sec)', 'fontsize',xfs);
    axis(AX)
end

x=.54;
y=.96;
w=.44;
h=.05;
annotation('textbox','String',{'Prior sampler'},...
    'FitHeightToText','on',...
    'LineStyle','none',...
    'fontsize',tfs,...
    'Position',[x y w h]); %[0.41 0.96 0.2398 0.05155]);

x=.41;
y=.47;
annotation('textbox','String',{'One observation ahead sampler'},...
    'FitHeightToText','on',...
    'LineStyle','none',...
    'fontsize',tfs,...
    'Position',[x y w h]); %[x,y, width, height]


% print to (color) eps
fig=figure(2);
wh=[7 7];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','~\Research\papers\BJ08\SimSampl')
save('SimSampl')

%% make a b&w fig
load('SimSampl')

n=1;
figure(3), cla, clf,

% gray    = [0.75 0.75 0.75];         %define gray
% col     = [1 0 0; .2 .2 1];          %define colors for mean
% ccol    = col+.4; ccol(ccol>1)=1;   %define colors for std
% 
% O           = R.F.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
% ONaNind     = find(~isfinite(O));
% Oind        = find(isfinite(O));
% O(ONaNind)  = [];
% finv        = ((P.k_d.*(P.beta-O))./(O-P.beta-P.alpha)).^(1/P.n);
% 
% nind    = find(R.n);                %find spike times
% nind    = nind(n);
% xmin    = Oind(find(Oind>nind,1)-1)-1;
% xmax    = Oind(find(Oind>nind,1))+1;
% xs      = Sim.tvec([xmin xmax]);%set the limits of the x-axis
% xind    = xmin:xmax;                %indices of x-axis
% ind     = Sim.T:-1:1;               %inverse index for 'fill' plots
% 
% %get min and max of calcium to normalize within plots
% 
% cshift=inf;
% for i=1:2
%     cmin(i) = min(min(R.C(xind)),min(min(S{i}.C(:,xind))));
%     cmax(i) = max(max(R.C(xind)),max(max(S{i}.C(:,xind))));
%     hmin(i) = min(min(P.omega*R.h(xind)),min(min(P.omega*S{i}.h(:,xind))));
%     hmax(i) = max(max(P.omega*R.h(xind)),max(max(P.omega*S{i}.h(:,xind))));
%     cshift  = min(min(min(S{i}.C(:,xind))),cshift);
%     M(i).hbar = sum(S{i}.w_b.*S{i}.h,1);
%     M(i).hvar = sum((repmat(M(i).hbar,Sim.N,1)-S{i}.h).^2)/Sim.N;
% end
% cmin    = min(cmin(:));
% cmax    = max(cmax(:));
% cdiff   = cmax-cmin;
% 
% hmin    = min(hmin(:));
% hmax    = max(hmax(:));
% hdiff   = hmax-hmin;
% 
% for i=1:2
%     C{i}=(S{i}.C(:,xind)-cmin)/cdiff;
%     hh{i}=(P.omega*S{i}.h(:,xind)-hmin)/hdiff;
% end
% 
% 
% %get forward means and variances
% for i=1:2
%     fCbar(i,:) = sum(S{i}.w_f.*S{i}.C,1);
%     fnbar(i,:) = sum(S{i}.w_f.*S{i}.n,1);
%     fnvar(i,:) = sum((repmat(fnbar(i,:),Sim.N,1)-S{i}.n).^2)/Sim.N;
% end
% 
% % set subfig sizes
% % w   = .41;      %width of subfigs
% % b   = .22;      %bottom of subfigs
% % h   = .3;       %height of subfigs
% % b2  = b+h+.05;  %bottom of 2nd row
% % l2  = .05+l1+w; %left of second col
% fs  = 12;       %default font size
% yfs = fs;       %ylabel font size
% xfs = fs;       %xlabel font size
% titfs = 12;   %title font size
% ticfs = fs;     %tick font size
% texfs = fs;     %text font size
% tfs   = 14;
% tl  = [.03 .03];%tick length
% sw  = 5;        %spike width
% bw  = .3;
% bw2 = .0015;
% xticks = xs(1):Sim.dt:xs(2);
% tx  = 1.435;
% ty  = 2.1;
% sp  = 1.2;
% lw  = 2;
% 
% Nrows=7;
% Ncols=2;
% AX  = [xs 0 1];
% interp = 'none';

% r2 = .5;
h  = .1;
w  = .35;
% l1 = .25;      %left of 1st col
% hs = .05;
% l2 = l1+w+hs;
% ms = 7;         %markersize
% ms2= 18;
% % yshift = .05;

% Sim.n = R.n;
% Sim.n(Sim.n==0) = NaN;
% obs = xind([2 5]);
for i=1:2

    % h particles
    subplot('Position',[l1 .34+(2-i)*r2 w h]) %[left bottom width height]    %     subplot(Nrows,Ncols,ii(1)),
    cla, hold on
    plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff,'color',gray,'LineWidth',2,'LineStyle','--')              %plot true calcium
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),hh{i}(nn,x)','.','Color','k','markersize',50*(S{i}.w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),hh{i}(:,1:length(xind))','Color','k')
    set(gca,'YTick',[0 1],'YTickLabel',[]); %,round([hmin hmax]*100)/100)
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    title('Particles','fontsize',titfs),
    ylab=ylabel([{'Weighted'}; {'Spike History'}],'fontsize',yfs,'Interpreter',interp);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

    % spike particles
    subplot('Position',[l1 .22+(2-i)*r2 w h]) %[left bottom width height]     subplot(Nrows,Ncols,ii(2)),
    cla, hold on
    for x=xind(1)+1:xind(end)
        bar(Sim.tvec(x),sum(S{i}.n(:,x))/5,'EdgeColor','k','FaceColor','k','BarWidth',bw2)                  %plot spike particles
        stem(Sim.tvec(x),Sim.n(x)-yshift,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
    end
    set(gca,'YTick',[0 1],'YTickLabel',[]); %[0 1])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    ylab=ylabel([{'Spike Train'}],'fontsize',yfs,'Interpreter',interp);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

    % calcium particles
    subplot('Position',[l1 .1+(2-i)*r2 w h]) %[left bottom width height]    % subplot(Nrows,Ncols,ii(3)), cla, hold on
    hold on
    plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2,'LineStyle','--')              %plot true calcium
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),C{i}(nn,x)','.','Color','k','markersize',50*(S{i}.w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),C{i}(:,1:length(xind))','Color','k')
    plot(Sim.tvec(obs),(R.C(obs)-cmin)/cdiff,'.k','Markersize',ms2)
    set(gca,'YTick',[0 1],'YTickLabel',[]); %,round((([cmin cmax])-cmin)*100)/100)
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
    xlabel('Time (sec)', 'fontsize',xfs);
    ylab=ylabel([{'Calcium'}; {'Concentration'}],'fontsize',yfs,'Interpreter',interp);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis(AX)
    box off


    % h dist
    subplot('Position',[l2 .34+(2-i)*r2 w h]) %[left bottom width height]%    subplot(Nrows,Ncols,ii(1)+1),
    cla, hold on
    ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.h);
    hfill=fill([Sim.tvec Sim.tvec(ind)],(P.omega*[ptiles(1,:) ptiles(2,ind)]-hmin)/hdiff,gray);
    set(hfill,'edgecolor',gray)
    plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff,'color',gray,'LineWidth',2,'LineStyle','--')              %plot true calcium
    plot(Sim.tvec,(P.omega*M(i).hbar-hmin)/hdiff,'linewidth',2,'color','k')
    set(gca,'YTick',[0 1],'YTickLabel',[])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    title('Inferred Distributions','fontsize',titfs),


    % spike plot
    subplot('Position',[l2 .22+(2-i)*r2 w h]) %[left bottom width height]    %subplot(Nrows,Ncols,ii(2)+1),
    cla, hold on
    ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.n);
    %     bar(Sim.tvec,R.n,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
    bar(Sim.tvec,ptiles(2,:),'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
    bar(Sim.tvec,M(i).nbar,'EdgeColor','k','FaceColor','k','BarWidth',bw)%plot true spikes
    stem(Sim.tvec,Sim.n-yshift,'Marker','v','MarkerSize',ms,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
    set(gca,'YTick',[0 1],'YTickLabel',[])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)


    % calcium plot
    subplot('Position',[l2 .1+(2-i)*r2 w h]) %[left bottom width height]    %subplot(Nrows,Ncols,ii(3)+1),
    cla, hold on
    ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.C);
    hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-cmin)/cdiff,gray);
    set(hfill,'edgecolor',gray)
    plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2,'LineStyle','--')
    plot(Sim.tvec,(M(i).Cbar-cmin)/cdiff,'linewidth',2,'color','k')
    plot(Sim.tvec(obs),(R.C(obs)-cmin)/cdiff,'.k','Markersize',ms2)
    set(gca,'YTick',[0 1],'YTickLabel',[])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
    xlabel('Time (sec)', 'fontsize',xfs);
    axis(AX)
end

x=.54;
y=.96;
w=.44;
h=.05;
annotation('textbox','String',{'Prior sampler'},...
    'FitHeightToText','on',...
    'LineStyle','none',...
    'fontsize',tfs,...
    'Position',[x y w h]); %[0.41 0.96 0.2398 0.05155]);

x=.41;
y=.47;
annotation('textbox','String',{'One observation ahead sampler'},...
    'FitHeightToText','on',...
    'LineStyle','none',...
    'fontsize',tfs,...
    'Position',[x y w h]); %[x,y, width, height]


% print to (color) eps
fig=figure(3);
wh=[7 7];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','~\Research\papers\BJ08\SimSampl_bw')

%% make old b&w fig
% 
% n=1;
% figure(3), cla, clf,
% 
% gray    = [0.75 0.75 0.75];         %define gray
% col     = [1 0 0; .2 .2 1];          %define colors for mean
% ccol    = col+.4; ccol(ccol>1)=1;   %define colors for std
% 
% O           = R.F.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);%let O be only observations at sample times
% ONaNind     = find(~isfinite(O));
% Oind        = find(isfinite(O));
% O(ONaNind)  = [];
% finv        = ((P.k_d.*(P.beta-O))./(O-P.beta-P.alpha)).^(1/P.n);
% 
% nind    = find(R.n);                %find spike times
% nind    = nind(n);
% xmin    = Oind(find(Oind>nind,1)-1)-1;
% xmax    = Oind(find(Oind>nind,1))+1;
% xs      = Sim.tvec([xmin xmax]);%set the limits of the x-axis
% xind    = xmin:xmax;                %indices of x-axis
% ind     = Sim.T:-1:1;               %inverse index for 'fill' plots
% 
% %get min and max of calcium to normalize within plots
% 
% cshift=inf;
% for i=1:2
%     cmin(i) = min(min(R.C(xind)),min(min(S{i}.C(:,xind))));
%     cmax(i) = max(max(R.C(xind)),max(max(S{i}.C(:,xind))));
%     hmin(i) = min(min(P.omega*R.h(xind)),min(min(P.omega*S{i}.h(:,xind))));
%     hmax(i) = max(max(P.omega*R.h(xind)),max(max(P.omega*S{i}.h(:,xind))));
%     cshift  = min(min(min(S{i}.C(:,xind))),cshift);
%     M(i).hbar = sum(S{i}.w_b.*S{i}.h,1);
%     M(i).hvar = sum((repmat(M(i).hbar,Sim.N,1)-S{i}.h).^2)/Sim.N;
% end
% cmin    = min(cmin(:));
% cmax    = max(cmax(:));
% cdiff   = cmax-cmin;
% 
% hmin    = min(hmin(:));
% hmax    = max(hmax(:));
% hdiff   = hmax-hmin;
% 
% for i=1:2
%     C{i}=(S{i}.C(:,xind)-cmin)/cdiff;
%     hh{i}=(P.omega*S{i}.h(:,xind)-hmin)/hdiff;
% end
% 
% 
% %get forward means and variances
% for i=1:2
%     fCbar(i,:) = sum(S{i}.w_f.*S{i}.C,1);
%     fnbar(i,:) = sum(S{i}.w_f.*S{i}.n,1);
%     fnvar(i,:) = sum((repmat(fnbar(i,:),Sim.N,1)-S{i}.n).^2)/Sim.N;
% end
% 
% % set subfig sizes
% % w   = .41;      %width of subfigs
% % b   = .22;      %bottom of subfigs
% % h   = .3;       %height of subfigs
% % b2  = b+h+.05;  %bottom of 2nd row
% % l2  = .05+l1+w; %left of second col
% fs  = 12;       %default font size
% yfs = fs;       %ylabel font size
% xfs = fs;       %xlabel font size
% titfs = 12;   %title font size
% ticfs = fs;     %tick font size
% texfs = fs;     %text font size
% tfs   = 14;
% tl  = [.03 .03];%tick length
% sw  = 5;        %spike width
% bw  = .3;
% bw2 = .0015;
% xticks = xs(1):Sim.dt:xs(2);
% tx  = 1.435;
% ty  = 2.1;
% sp  = 1.2;
% lw  = 2;
% 
% Nrows=7;
% Ncols=2;
% AX  = [xs 0 1];
% interp = 'none';
% 
% r2 = .5;
% h  = .1;
% l1 = .25;      %left of 1st col
% hs = .05;
% w  = .35;
% l2 = l1+w+hs;
% 
% for i=1:2
% 
%     %     if i==1
%     %         ii=[1 3 5];
%     %     else
%     %         ii=[9 11 13];
%     %     end
%     %
%     % h particles
%     subplot('Position',[l1 .34+(2-i)*r2 w h]) %[left bottom width height]
%     %     subplot(Nrows,Ncols,ii(1)),
%     cla, hold on
%     plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff,'color',gray,'LineWidth',2)              %plot true calcium
%     for x=1:length(xind)
%         for nn=1:Sim.N
%             plot(Sim.tvec(xind(x)),hh{i}(nn,x)','.','Color','k','markersize',50*(S{i}.w_f(nn,x)))                   %plot calcium particles
%         end
%     end
%     plot(Sim.tvec(xind),hh{i}(:,1:length(xind))','Color','k')
%     set(gca,'YTick',[0 1],'YTickLabel',[]); %,round([hmin hmax]*100)/100)
%     set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
%     axis(AX)
%     %     if i==1,                                                                        %bottom left
%     title('Particles','fontsize',titfs),
%     %     end
%     ylab=ylabel([{'Weighted'}; {'Spike History'}],'fontsize',yfs,'Interpreter',interp);
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% 
%     % spike particles
%     subplot('Position',[l1 .22+(2-i)*r2 w h]) %[left bottom width height]     subplot(Nrows,Ncols,ii(2)),
%     cla, hold on
%     for x=xind(1)+1:xind(end)
%         bar(Sim.tvec(x),R.n(x),'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw2)%plot true spikes
%         bar(Sim.tvec(x),sum(S{i}.n(:,x))/5,'EdgeColor','k','FaceColor','k','BarWidth',bw2)                  %plot spike particles
%     end
%     set(gca,'YTick',[0 1],'YTickLabel',[]); %[0 1])
%     set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
%     axis(AX)
%     ylab=ylabel([{'Spike Train'}],'fontsize',yfs,'Interpreter',interp);
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% 
%     % calcium particles
%     subplot('Position',[l1 .1+(2-i)*r2 w h]) %[left bottom width height]
%     % subplot(Nrows,Ncols,ii(3)), cla, hold on
%     plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2)              %plot true calcium
%     for x=1:length(xind)
%         for nn=1:Sim.N
%             plot(Sim.tvec(xind(x)),C{i}(nn,x)','.','Color','k','markersize',50*(S{i}.w_f(nn,x)))                   %plot calcium particles
%         end
%     end
%     plot(Sim.tvec(xind),C{i}(:,1:length(xind))','Color','k')
%     set(gca,'YTick',[0 1],'YTickLabel',[]); %,round((([cmin cmax])-cmin)*100)/100)
%     set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
%     xlabel('Time (sec)', 'fontsize',xfs);
%     ylab=ylabel([{'Calcium'}; {'Concentration'}],'fontsize',yfs,'Interpreter',interp);
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
%     axis(AX)
%     box off
% 
% 
%     % h dist
%     subplot('Position',[l2 .34+(2-i)*r2 w h]) %[left bottom width height]
%     %    subplot(Nrows,Ncols,ii(1)+1),
%     cla, hold on
%     ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.h);
%     hfill=fill([Sim.tvec Sim.tvec(ind)],(P.omega*[ptiles(1,:) ptiles(2,ind)]-hmin)/hdiff,gray);
%     set(hfill,'edgecolor',gray)
%     plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff,'color',gray,'LineWidth',2)              %plot true calcium
%     plot(Sim.tvec,(P.omega*M(i).hbar-hmin)/hdiff,'linewidth',2,'color','k')
%     set(gca,'YTick',[0 1],'YTickLabel',[])
%     set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
%     axis(AX)
%     %     if i==1,                                                                        %bottom left
%     title('Inferred Distributions','fontsize',titfs),
%     %     end
% 
%     % spike plot
%     subplot('Position',[l2 .22+(2-i)*r2 w h]) %[left bottom width height]
%     %subplot(Nrows,Ncols,ii(2)+1),
%     cla, hold on
%     ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.n);
%     bar(Sim.tvec,R.n,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
%     bar(Sim.tvec,ptiles(2,:),'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
%     bar(Sim.tvec,M(i).nbar,'EdgeColor','k','FaceColor','k','BarWidth',bw)%plot true spikes
%     set(gca,'YTick',[0 1],'YTickLabel',[])
%     set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
%     axis(AX)
% 
% 
%     % calcium plot
%     subplot('Position',[l2 .1+(2-i)*r2 w h]) %[left bottom width height]
%     %subplot(Nrows,Ncols,ii(3)+1),
%     cla, hold on
%     ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.C);
%     hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-cmin)/cdiff,gray);
%     set(hfill,'edgecolor',gray)
%     plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2)
%     plot(Sim.tvec,(M(i).Cbar-cmin)/cdiff,'linewidth',2,'color','k')
%     set(gca,'YTick',[0 1],'YTickLabel',[])
%     set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
%     xlabel('Time (sec)', 'fontsize',xfs);
%     axis(AX)
% end
% 
% x=.54;
% y=.96;
% w=.44;
% h=.05;
% annotation('textbox','String',{'Prior sampler'},...
%     'FitHeightToText','on',...
%     'LineStyle','none',...
%     'fontsize',tfs,...
%     'Position',[x y w h]); %[0.41 0.96 0.2398 0.05155]);
% 
% x=.41;
% y=.47;
% annotation('textbox','String',{'One observation ahead sampler'},...
%     'FitHeightToText','on',...
%     'LineStyle','none',...
%     'fontsize',tfs,...
%     'Position',[x y w h]); %[x,y, width, height]
% 
% 
% % print to (color) eps
% fig=figure(3);
% wh=[7 7];
% set(fig,'PaperPosition',[0 11-wh(2) wh]);
% print('-depsc','~\Research\papers\BJ08\SimSampl_bw')