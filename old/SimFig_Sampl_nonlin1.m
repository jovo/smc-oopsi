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
Sim.freq    = 5;                        %frequency of observations
Sim.Nsec    = 2.5;                      %# of sec
Sim.StimDim = 1;                        %# of stimulus dimensions
Sim.M       = 0;                        %number of spike history terms
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

% get "real" data
R.n         = zeros(1,Sim.T);   %spike times
R.C         = P.C_init*ones(1,Sim.T);   %initialize calcium
epsilon_c   = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);%generate noise on calcium
spt         = [3*Sim.freq:Sim.freq:50 81 131 181 231 281 331 381 431 481];    %forced spike times
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
n=11;
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
shift   = 20;
xs      = [Sim.tvec(nind(1)-shift) Sim.tvec(nind(1)+shift)];%set the limits of the x-axis
xmin    = find(Sim.tvec>xs(1),1);   %index of lower limit on x
xmax    = find(Sim.tvec>xs(2),1)-1; %index of upper limit on x
xind    = xmin:xmax;                %indices of x-axis
ind     = Sim.T:-1:1;               %inverse index for 'fill' plots

%get min and max of calcium to normalize within plots
for i=1:2
    cmin(i) = min(min(M(i).Cbar(xind)-sqrt(M(i).Cvar(xind))));
    cmax(i) = max(max(M(i).Cbar(xind)+sqrt(M(i).Cvar(xind))));
end
cmin    = min(cmin(:));
cmax    = max(cmax(:))-P.C_0;
AX      = [xs 0 2*cmax];                 %axes

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
bw  = 1;
xticks = AX(1):Sim.dt*Sim.freq*2:AX(2);
tx  = 1.435;
ty  = 2.1;

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
        %         text(tx,ty,'(A)','fontsize',texfs);
    end
    if i==2,                                                                        %top right plot
        subplot('Position',[l1 b w h]), %subplot(2,2,3)
        cla, hold on
        ylab=ylabel({'Conditional'; 'Sampler'});
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',xticks-AX(1),'TickLength',tl,'fontsize',ticfs,'YTick',[])        %         text(tx,ty,'(C)','fontsize',texfs);
    end
    set(ylab,'fontsize',yfs)%'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    axis(AX)
    set(gca,'YTickLabel',[])

    % calcium particles
    plot(Sim.tvec,S{i}.C'-P.C_0+cmax,'Color',ccol(i,:))                   %plot calcium particles
    plot(Sim.tvec,R.C-P.C_0+cmax,'color',gray,'LineWidth',2)              %plot true calcium
    plot(Oind*Sim.dt,finv+cmax,'.k','LineWidth',1,'markersize',7)      %plot observations
    plot(Sim.tvec,cmax*ones(size(Sim.tvec)),'k')

    % spike particles
%     stem(Sim.tvec,R.n*cmax,'Marker','none','Color',gray,'LineWidth',sw)                  %plot true spikes
	bar(Sim.tvec,R.n*cmax,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
%  	bar(Sim.tvec,S{i}.n'*cmax,'EdgeColor',ccol(i,:),'FaceColor',ccol(i,:),'BarWidth',.1)%plot true spikes
    stem(Sim.tvec,S{i}.n'*cmax,'Marker','none','Color',ccol(i,:),'LineWidth',1)                  %plot true spikes
    plot(Sim.tvec,zeros(size(Sim.tvec)),'k')

    %% plot inferred distributions
    if i==1,                                                                        %bottom left
        subplot('Position',[l2 b2 w h]), %subplot(2,2,2),
        cla, hold on
        title('Inferred Distributions','fontsize',titfs),
        set(gca,'XTickLabel',[])
        set(gca,'XTick',xticks,'XTickLabel',[],'TickLength',tl,'YTick',[])
        %         text(tx,ty,'(B)','fontsize',texfs);
    else                                                                            %bottom right
        subplot('Position',[l2 b w h]),%subplot(2,2,4),
        cla, hold on
        xlabel('Time (sec)', 'fontsize',xfs);
        set(gca,'XTick',xticks,'XTickLabel',xticks-AX(1),'TickLength',tl,'fontsize',ticfs,'YTick',[])
        %         text(tx,ty,'(D)','fontsize',texfs);
    end
    axis(AX)
    set(gca,'YTickLabel',[])

    % calcium plot
    %hfill=fill([Sim.tvec Sim.tvec(ind)],[M(i).Cbar-sqrt(M(i).Cvar) M(i).Cbar(ind)+sqrt(M(i).Cvar(ind))]-P.C_0+cmax,ccol(i,:));
    ptiles = GetPercentiles([.25 .75],S{i}.w_b,S{i}.C);
    hfill=fill([Sim.tvec Sim.tvec(ind)],[ptiles(1,:) ptiles(2,ind)]+cmax-P.C_0,ccol(i,:));
    set(hfill,'edgecolor',ccol(i,:))
    plot(Sim.tvec,R.C-P.C_0+cmax,'color',gray,'LineWidth',2)
    plot(Sim.tvec,M(i).Cbar-P.C_0+cmax,'linewidth',2,'color',col(i,:))
    plot(Oind*Sim.dt,finv+cmax,'.k','LineWidth',1,'markersize',7)      %plot observations
    plot(Sim.tvec,cmax*ones(size(Sim.tvec)),'k')

    % spike plot
    BarVar=M(i).nbar+M(i).nvar;                                                     %make var of spikes not exceed 1
    BarVar(BarVar>1)=1;
	bar(Sim.tvec,R.n*cmax,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
	bar(Sim.tvec,BarVar*cmax,'EdgeColor',ccol(i,:),'FaceColor',ccol(i,:),'BarWidth',bw)%plot true spikes
	bar(Sim.tvec,M(i).nbar*cmax,'EdgeColor',col(i,:),'FaceColor',col(i,:),'BarWidth',bw)%plot true spikes
%     stem(Sim.tvec,BarVar*cmax,'Marker','none','Color',ccol(i,:),'LineWidth',sw)          %plot forward spike var
%     stem(Sim.tvec,M(i).nbar*cmax,'Marker','none','Color',col(i,:),'LineWidth',sw)        %plot forward spike mean
end


% print to (color) eps
fig=figure(2);
wh=[7 4.5];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\sampl2