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
    Sim.pf=i-1;
    [S(i) M(i)] = smc_em_bern_FoBaMo_v5(Sim,R,P);
end
%% make a fig
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

for i=1:2
    cshift=inf;
    cmin(i) = min(min(R.C(xind)),min(min(S(i).C(:,xind))));
    cmax(i) = max(max(R.C(xind)),max(max(S(i).C(:,xind))));
    cshift(i)  = min(min(min(S(i).C(:,xind))),cshift);

    hmin(i) = min(min(P.omega*R.h(xind)),min(min(P.omega*S(i).h(:,xind))));
    hmax(i) = max(max(P.omega*R.h(xind)),max(max(P.omega*S(i).h(:,xind))));
    M(i).hbar = sum(S(i).w_b.*S(i).h,1);
    M(i).hvar = sum((repmat(M(i).hbar,Sim.N,1)-S(i).h).^2)/Sim.N;
end


cmin    = min(cmin(:));
cmax    = max(cmax(:));
cdiff   = cmax-cmin;
hmin    = min(hmin(:));
hmax    = max(hmax(:));
hdiff   = hmax-hmin;

for i=1:2
    C(i)=(S(i).C(:,xind)-cmin)/cdiff;
    hh=(P.omega*S(i).h(:,xind)-hmin)/hdiff;
end


%get forward means and variances
fCbar = sum(S(i).w_f.*S(i).C,1);
fnbar = sum(S(i).w_f.*S(i).n,1);
fnvar = sum((repmat(fnbar(i,:),Sim.N,1)-S(i).n).^2)/Sim.N;

% set subfig sizes
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

Nrows=3;
Ncols=2;
AX  = [xs 0 1];
ii=[1 3 5];

for i=1:2
    % h particles
    subplot(Nrows,Ncols,ii(1)+i-1), cla, hold on
    plot(Sim.tvec(xind),(P.omega*R.h(xind)-hmin)/hdiff,'color',gray,'LineWidth',2)              %plot true calcium
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),hh(nn,x)','.','Color',col(2,:),'markersize',50*(S(i).w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),hh(:,1:length(xind))','Color',col(2,:))
    set(gca,'YTick',[0 1],'YTickLabel',round([hmin hmax]*100)/100)
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    title('Particles','fontsize',titfs,'Interpreter','latex'),
    ylabel({'$\omega h_t$ (a.u.)'},'fontsize',yfs,'Interpreter','latex')


    % spike particles
    subplot(Nrows,Ncols,ii(2)+i-1), cla, hold on
    for x=xind(1)+1:xind(end)
        bar(Sim.tvec(x),R.n(x),'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw2)%plot true spikes
        bar(Sim.tvec(x),sum(S(i).n(:,x))/5,'EdgeColor',col(2,:),'FaceColor',col(2,:),'BarWidth',bw2)                  %plot spike particles
    end
    set(gca,'YTick',[0 1],'YTickLabel',[0 1])
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
    axis(AX)
    ylabel({'$n_t$ $(\#)$'},'fontsize',yfs,'Interpreter','latex')

    % calcium particles
    subplot(Nrows,Ncols,ii(3)+i-1), cla, hold on
    plot(Sim.tvec(xind),(R.C(xind)-cmin)/cdiff,'color',gray,'LineWidth',2)              %plot true calcium
    plot(Sim.tvec(xind(2)), (R.C(xind(2))-cmin)/cdiff,'ok','LineWidth',2,'markersize',11)
    plot(Sim.tvec(xind(end)-1), (R.C(xind(end)-1)-cmin)/cdiff,'ok','LineWidth',2,'markersize',11)
    for x=1:length(xind)
        for nn=1:Sim.N
            plot(Sim.tvec(xind(x)),C(i)(nn,x)','.','Color',col(2,:),'markersize',50*(S(i).w_f(nn,x)))                   %plot calcium particles
        end
    end
    plot(Sim.tvec(xind),C(i)(:,1:length(xind))','Color',col(2,:))
    set(gca,'YTick',[0 1],'YTickLabel',round((([cmin cmax])-cmin)*10)/10)
    set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
    xlabel('Time (sec)', 'fontsize',xfs,'Interpreter','latex');
    ylabel([{'[Ca$^{2+}]_t$ ($\mu$M)'}],'Interpreter','latex','fontsize',yfs)
    axis(AX)
end

% % h dist
% subplot(Nrows,Ncols,ii(1)+1), cla, hold on
% ptiles = GetPercentiles([.25 .75],S(i).w_b,S(i).h);
% hfill=fill([Sim.tvec Sim.tvec(ind)],(P.omega*[ptiles(1,:) ptiles(2,ind)]-hmin)/hdiff,ccol(i,:));
% set(hfill,'edgecolor',ccol(i,:))
% plot(Sim.tvec,(P.omega*R.h-hmin)/hdiff,'color',gray,'LineWidth',2)              %plot true calcium
% plot(Sim.tvec,(P.omega*M(i).hbar-hmin)/hdiff,'linewidth',2,'color',col(2,:))
% set(gca,'YTick',[0 1],'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
% axis(AX)
% if i==1,                                                                        %bottom left
%     title('Inferred Distributions','fontsize',titfs),
% end
%
% % spike plot
% subplot(Nrows,Ncols,ii(2)+1), cla, hold on
% ptiles = GetPercentiles([.25 .75],S(i).w_b,S(i).n);
% %     BarVar=M(i).nbar+M(i).nvar; BarVar(BarVar>1)=1;
% bar(Sim.tvec,R.n,'EdgeColor',gray,'FaceColor',gray,'BarWidth',bw)%plot true spikes
% bar(Sim.tvec,ptiles(2,:),'EdgeColor',ccol(i,:),'FaceColor',ccol(i,:),'BarWidth',bw)%plot true spikes
% bar(Sim.tvec,M(i).nbar,'EdgeColor',col(2,:),'FaceColor',col(2,:),'BarWidth',bw)%plot true spikes
% set(gca,'YTick',[0 1],'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[])
% axis(AX)
%
%
% % calcium plot
% subplot(Nrows,Ncols,ii(3)+1), cla, hold on
% ptiles = GetPercentiles([.25 .75],S(i).w_b,S(i).C);
% hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-cmin)/cdiff,ccol(i,:));
% set(hfill,'edgecolor',ccol(i,:))
% plot(Sim.tvec,(R.C-cmin)/cdiff,'color',gray,'LineWidth',2)
% plot(Sim.tvec,(M(i).Cbar-cmin)/cdiff,'linewidth',2,'color',col(2,:))
% %     plot(Oind*Sim.dt,(finv-cmin)/cdiff,'.k','LineWidth',1,'markersize',7)      %plot observations
% set(gca,'YTick',[0 1],'YTickLabel',[])
% set(gca,'XTick',Sim.tvec(xind),'XTickLabel',[{''};{'u'}; {''}; {''}; {'v'}])
% xlabel('Time (sec)', 'fontsize',xfs);
% axis(AX)

% annotation('textbox','String',{'Prior Sampler'},...
%     'FitHeightToText','on',...
%     'LineStyle','none',...
%     'fontsize',tfs,...
%     'Position',[0.38 0.45 0.2398 0.05155]);
%
% annotation('textbox','String',{'Prior Sampler'},...
%     'FitHeightToText','on',...
%     'LineStyle','none',...
%     'fontsize',tfs,...
%     'Position',[0.41 0.96 0.2398 0.05155]);

% print to (color) eps
fig=figure(2);
wh=[5 7];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\working_copies\neur_ca_imag\trunk\columbia_talk\sampl_comp