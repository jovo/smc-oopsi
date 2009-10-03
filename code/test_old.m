clc, clear V P
%% 2) set simulation metadata

% V.T_o     = 500;                                    % # of time steps
V.dt      = 1/30;                            % time step size
% V.freq    = 1;                                    % # of time steps between observations
% V.T       = V.T_o*V.freq;                     % # of observations
% V.Nsec    = V.T*V.dt;                         % # of actual seconds
% V.tvec    = V.dt:V.dt:V.Nsec;               % vector of times
% V.N       = 100;                                  % # of particles
% V.M       = 0;                                    % number of spike history dimensions
% V.pf      = 1;                                    % use conditional sampler (not prior) when possible
% V.StimDim = 1;                                    % # of stimulus dimensions
% V.x       = ones(1,V.T);                        % stimulus
% V.smc_plot = 1;

% V.Mstep   = 1;                                 % do M-step
% V.est_c = 1;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
% V.est_n = 1;                                 % whether to estimate rate governing parameters {b,k}
% V.est_h = 0;                                % whether to estimate spike history parameters {h}
% V.est_F = 1;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
% V.smc_iter_max = 2;                                    % max # of EM iterartions
% V.fast = 0;
% V.Scan    = 1;

%% 3) initialize parameters

% initialize particle filter parameters
P.k         = -0.61; %log(-log(1-sum(D{j}.n)/V.T)/V.dt);% linear filter
P.tau_c     = 1.09;%V.dt/x(2);
P.A         = 51.4;
P.C_0       = 1.99;                                   % baseline [Ca++]
P.C_init    = 0.02;                                 % initial [Ca++]
P.sigma_c   = 6.25;
P.n         = 1;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = .87;                                    % F_max
P.beta      = 0.09;                                    % F_min
P.gamma     = 0;                                 % scaled variance
P.zeta      = 0.01;                                 % constant variance

%% 4) infer spikes and estimate parameters

[smc.E smc.P smc.V] = GOOPSI_main_v1_0(F,V,P);

%% 5) plot results

fig=figure(3); clf,
V.T     = length(F);
% if V.fast_do==1, V.dt=fast.V.dt; else  V.dt=smc.V.dt; end
% if (V.fast_do==1 && V.smc_do==1), nrows=3; else nrows=2; end
nrows   = 2;
gray    = [.75 .75 .75];            % define gray color
col   = [1 0 0; 0.2 0.2 1];         % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter   = 'tex';                    % interpreter for axis labels
xlims   = [45 V.T-2];               % xmin and xmax for current plot
fs      = 12;                       % font size
ms      = 5;                        % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
xticks  = 0:1/V.dt:V.T;             % XTick positions
skip    = round(length(xticks)/5);
xticks  = xticks(1:skip:end);
tvec_o=xlims(1):xlims(2);
i=0;
if isfield(V,'true_n'), spt = find(V.true_n==1); end               % find spike times

% plot fluorescence data
i=i+1; h(i)=subplot(nrows,1,i); hold on
plot(tvec_o,z1(F(tvec_o)),'.-k','LineWidth',2,'MarkerSize',ms*.75);
% stem(V.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])


% % plot fast-oopsi output
% if V.fast_do==1
%     i=i+1; h(i)=subplot(nrows,1,i); hold on,
%     n_fast=fast.n/max(fast.n);
%     spts=find(n_fast>1e-3);
%     stem(spts,n_fast(spts),'Marker','none','LineWidth',sw,'Color','k')
%     if isfield(V,'true_n'),
%         stem(spt,V.true_n(spt)+0.1,'Marker','v','MarkerSize',ms,'LineStyle','none','Color',gray,'MarkerFaceColor',gray,'MarkerEdgeColor',gray)
%     end
%     axis([xlims 0 1.1])
%     hold off,
%     ylab=ylabel([{'Fast'}; {'Filter'}],'Interpreter',inter,'FontSize',fs);
%     set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
%     set(gca,'YTick',0:2,'YTickLabel',[])
%     set(gca,'XTick',xticks,'XTickLabel',[])
%     box off
% end
% 
% % plot smc-oopsi output
% if V.smc_do == 1
    i=i+1; h(i)=subplot(nrows,1,i); hold on,
    spts=find(smc.E.nbar>1e-3);
    stem(spts,smc.E.nbar(spts),'Marker','none','LineWidth',sw,'Color','k')
    if isfield(V,'true_n'),
        stem(spt,V.true_n(spt)+0.1,'Marker','v','MarkerSize',ms,'LineStyle','none','Color',gray,'MarkerFaceColor',gray,'MarkerEdgeColor',gray)
    end
    axis([xlims 0 1.1])
    hold off,
    ylab=ylabel([{'Particle'}; {'Filter'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    box off
% end

% label last subplot
set(gca,'XTick',xticks,'XTickLabel',round(xticks*V.dt*100)/100)
xlabel('Time (sec)','FontSize',fs)
linkaxes(h,'x')

% print fig
wh=[7 3];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
% print('-depsc',V.name_fig)
% print('-dpdf',V.name_fig)
% saveas(fig,V.name_fig)
