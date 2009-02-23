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
Sim.M       = 0;                        %number of spike history terms
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
P.gamma     = 1e-5;                     %var gainh
P.zeta      = .1e-5;                    %var offset
P.a         = Sim.dt/P.tau_c;
P.A         = .05;
P.sigma_c   = .1;

% get "real" data
R.n         = zeros(1,Sim.T);   %spike times
R.C         = P.C_init*ones(1,Sim.T);   %initialize calcium
epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
spt         = [3*Sim.freq:Sim.freq:50 81 131 181 231 281 331 381 431 481]-1;    %forced spike times
spt         = spt(spt<Sim.T);
R.n(spt)    = 1;                %force spikes

for t=2:Sim.T                   %update calcium
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
n=16;
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
    cshift  = min(min(min(S{i}.C(:,xind))),cshift);
end
cmin    = min(cmin(:));
cmax    = max(cmax(:));
cdiff   = cmax-cmin;
AX      = [xs 0 2*(cmax-cmin)];                 %axes

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

% [foo ind]=find(R.n==0);
% if ~isempty(ind), R.n(ind)=NaN; end

% make plots
for i=[2 1]

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
    axis(AX)

    % calcium particles
    plot(Sim.tvec,R.C+cdiff-cshift,'color',gray,'LineWidth',2)              %plot true calcium
    for x=xind
        for nn=1:Sim.N
            plot(Sim.tvec(x),S{i}.C(nn,x)'+cdiff-cshift,'.','Color',col(i,:),'markersize',10*exp(S{i}.w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),S{i}.C(:,xind)'+cdiff-cshift,'Color',col(i,:))
    plot(Oind*Sim.dt,finv+cdiff-cshift,'.k','LineWidth',1,'markersize',7)      %plot observations
    plot(Sim.tvec,cdiff*ones(size(Sim.tvec)),'k')

    % spike particles
    for x=xind
    	bar(Sim.tvec(x),R.n(x)*cdiff,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw2)%plot true spikes
        bar(Sim.tvec(x),cdiff*sum(S{i}.n(:,x))/5,'EdgeColor',col(i,:),'FaceColor',col(i,:),'BarWidth',bw2)                  %plot spike particles
    end
    plot(Sim.tvec,zeros(size(Sim.tvec)),'k')


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