

%let O be only observations at sample times
O           = R.F.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
ONaNind     = find(~isfinite(O));
Oind        = find(isfinite(O));
O(ONaNind)  = [];

% do some x-axis stuff
xmin    = Oind(2);
xmax    = Oind(end-6);
xind    = xmin:xmax;
xs      = [Sim.tvec(xmin) Sim.tvec(xmax)]+Sim.dt;
spikeAX = [xs 0 1];
hidAX   = [xs 0 2];
xticks  = Sim.tvec(Oind(1))+Sim.dt:Sim.dt*Sim.freq*2:Sim.tvec(Oind(end));
xticks  = xticks-xticks(1);

% do some other stuff
M.pbar = sum(S.w_b.*S.p,1);
M.pvar = sum((repmat(M.pbar,Sim.N,1)-S.p).^2)/Sim.N;
pmin    = min(min(min(R.p(xind),min(M.pbar(xind)-sqrt(M.pvar(xind))))));
pmax    = max(max(max(R.p(xind),max(M.pbar(xind)+sqrt(M.pvar(xind))))));

M.hbar = sum(S.w_b.*S.h,1);
M.hvar = sum((repmat(M.hbar,Sim.N,1)-S.h).^2)/Sim.N;
hmin    = min(min(min(P.omega*R.h(xind),P.omega*min(M.hbar(xind)-sqrt(M.hvar(xind))))));
hmax    = max(max(max(P.omega*R.h(xind),P.omega*max(M.hbar(xind)+sqrt(M.hvar(xind))))));

cmin    = min(min(min(R.C(xind),min(M.Cbar(xind)-sqrt(M.Cvar(xind))))));
cmax    = max(max(max(R.C(xind),max(M.Cbar(xind)+sqrt(M.Cvar(xind))))));

gray    = [0.75 0.75 0.75];
col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8; ccol(ccol>1)=1;
ind     = Sim.T:-1:1;
sw      = 4;
lw      = 2;
fs      = 14;
tfs     = 14;

figure(1), clf, 
Nrows=8;
Ncols=2;
ii=7; 
set(gcf, 'color', 'w');

% external stimulus
subplot('Position',[0.18 0.838 0.285 0.08]), cla, hold on
imagesc(Sim.x)
colormap('gray')
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xmin xmax 0.5 5.5])
% ylab=ylabel({'External'; 'Stimulus'});
ylab=ylabel('$x_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% title('Input','fontsize',tfs)

% filtered stimulus
subplot('Position',[0.18 0.73 0.285 0.08]), cla, hold on, box on
plot(Sim.tvec,P.k'*Sim.x,'color','k','linewidth',lw)
% set(gca,'XTick',xticks,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% set(gca,'YTick',[])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(P.k'*Sim.x) max(P.k'*Sim.x)])
% ylab=ylabel({'Linearly Filtered'; 'Stimulus'});
ylab=ylabel('$k^T x_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% xlab=xlabel('Time (sec)');

text(0.45, -28,'Observation State','fontsize',tfs)
text(0.58, 8,'Input','fontsize',tfs)

% P(spike) 
i=ii; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,(R.p-pmin)/(pmax-pmin),'Color',gray,'LineWidth',lw)
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel('$p_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike train
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel('$n_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike history
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,(P.omega*R.h-hmin)/(hmax-hmin),'Color',gray,'LineWidth',lw)%(R.h-hmin)/(hmax-hmin)
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel('$h_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% calcium
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray,'LineWidth',lw)%(R.C-cmin)/(cmax-cmin)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs 0 1])
% set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
ylab=ylabel('[Ca$^{2+}]_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% xlab=xlabel('Time (sec)');

% observation state
subplot('Position',[0.18 0.05 0.285 0.08]), cla, hold on, box on
plot(Oind*Sim.dt,O,'.k','LineWidth',1,'markersize',7)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(O(2:end-6)) max(O(2:end-6))])
xlabel('Fluorescence Observations');
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% title('Observation State','fontsize',tfs)
ylab=ylabel('$F_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlab=xlabel('Time (sec)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% external stimulus
subplot('Position',[0.62 0.838 0.285 0.08]), cla, hold on
imagesc(Sim.x)
colormap('gray')
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xmin xmax 0.5 5.5])
ylab=ylabel({'External      '; 'Stimulus     '});
% ylab=ylabel('$x_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% title('Input','fontsize',tfs)

% filtered stimulus
subplot('Position',[0.62 0.73 0.285 0.08]), cla, hold on, box on
plot(Sim.tvec,P.k'*Sim.x,'color','k','linewidth',lw)
% set(gca,'XTick',xticks,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% set(gca,'YTick',[])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(P.k'*Sim.x) max(P.k'*Sim.x)])
ylab=ylabel({'Linearly Filtered'; 'Stimulus    '});
% ylab=ylabel('$k^T x_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% xlab=xlabel('Time (sec)');

% P(spike)
i=ii+1; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,S.p);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-pmin)/(pmax-pmin),ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(R.p-pmin)/(pmax-pmin),'Color',gray,'LineWidth',lw)
plot(Sim.tvec,(M.pbar-pmin)/(pmax-pmin),'linewidth',2,'color',col(2,:))
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel({'Probability   ';'of spiking    '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike train
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
BarVar=M.nbar+M.nvar;
BarVar(BarVar>1)=1;
stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',sw)
stem(Sim.tvec,M.nbar,'Marker','none','Color',col(2,:),'LineWidth',sw)
axis([xs min(R.n) max(R.n)])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel({'Spike        '; 'Train         '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike history
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,P.omega*S.h);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-hmin)/(hmax-hmin),ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(P.omega*R.h-hmin)/(hmax-hmin),'Color',gray,'LineWidth',lw)%(R.h-hmin)/(hmax-hmin)
plot(Sim.tvec,(P.omega*M.hbar-hmin)/(hmax-hmin),'linewidth',2,'color',col(2,:))
axis([xs 0 1])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
ylab=ylabel({'Spike History ';'Term        '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% calcium
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,S.C);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-cmin)/(cmax-cmin),ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'Color',gray,'LineWidth',lw)%(R.C-cmin)/(cmax-cmin)
plot(Sim.tvec,(M.Cbar-cmin)/(cmax-cmin),'linewidth',2,'color',col(2,:))
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs 0 1])
ylab=ylabel({'Calcium     '; 'Concentration '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% xlab=xlabel('Time (sec)');

text(-0.60, 5.6,'True Hidden States','fontsize',tfs)
text(+0.06, 5.6,'Inferred Hidden States','fontsize',tfs)

subplot('Position',[0.62 0.05 0.285 0.08]), cla, hold on, box on
plot(Oind*Sim.dt,O,'.k','LineWidth',1,'markersize',7)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(O(2:end-6)) max(O(2:end-6))])
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% title('Observation State','fontsize',tfs)
ylab=ylabel({'Fluorescence  ';'Observations  '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlab=xlabel('Time (sec)');


% subplot('Position',[15 0.05 0.01 0.1]), cla, hold on, box off
% axes('ZColor',[1 1 1],'YColor',[1 1 1],'XColor',[1 1 1],'Position',[0.786 0.05 0.1 0.1]);
% ylab=ylabel({'Fluorescence'; 'Observations'});
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color','k')
% ['Inferred',sprintf('\n'),'Hidden States']
% text(+0.18, -1.9,['Fluorescence',sprintf('\n'),'Observations'])
% text(+0.18, 8.3,['Linearly Filtered',sprintf('\n'),'Stimulus'])
% text(+0.18, 10.4,['External',sprintf('\n'),'Stimulus'])

% make eps
fig=figure(1);
wh=[7 7];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\schem3