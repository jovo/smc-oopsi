

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
sw      = 1.8;
lw      = 1.8;
fs      = 14;
tfs     = 14;

figure(1), clf, 
Nrows=8;
Ncols=2;
ii=7; 
set(gcf, 'color', 'w');

% external stimulus
subplot('Position',[0.166 0.838 0.285 0.08]), cla, hold on
imagesc(Sim.x)
colormap('gray')
set(gca,'YTick',[1 Sim.StimDim],'YTickLabel',[1 Sim.StimDim],'XTick',[],'FontName','Times New Roman');
axis([xmin xmax 0.5 5.5])
% ylab=ylabel({'External'; 'Stimulus'});
ylab=ylabel([{'\quad \quad             $x_t$'}; {'(dimensions)'}],'Interpreter','latex','VerticalAlignment','middle','FontName','Times New Roman');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% title('Input','fontsize',tfs)

% filtered stimulus
subplot('Position',[0.166 0.73 0.285 0.08]), cla, hold on, box on
plot(Sim.tvec,P.k'*Sim.x,'color','k','linewidth',lw)
% set(gca,'XTick',xticks,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% set(gca,'YTick',[])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(P.k'*Sim.x) max(P.k'*Sim.x)])
% ylab=ylabel({'Linearly Filtered'; 'Stimulus'});
ylab=ylabel([{'$k^T x_t$'}; {'(a.u.)'}],'Interpreter','latex','VerticalAlignment','middle','FontName','Times New Roman');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% xlab=xlabel('Time (sec)');
set(gca,'YTick',[-1.5 1.5],'YTickLabel',{'-1.5','+1.5'},'FontName','Times New Roman');

text(0.1, -6.2,'True Hidden States','fontsize',tfs,'FontName','Times New Roman');
text(+0.76, -6.2,'Inferred Hidden States','fontsize',tfs,'FontName','Times New Roman');
text(0.48, -28,'Observation State','fontsize',tfs,'FontName','Times New Roman');
text(0.58, 8.2,'Input','fontsize',tfs,'FontName','Times New Roman');

% P(spike) 
i=ii; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,(R.p-pmin)/(pmax-pmin),'Color',gray,'LineWidth',lw)
axis([xs 0 1])
set(gca,'XTickLabel',[],'XTick',xticks,'FontName','Times New Roman');
set(gca,'YTick',[0 1],'FontName','Times New Roman');
ylab=ylabel([{'  $p_t$'}; {'($\%$)'}],'Interpreter','latex','VerticalAlignment','middle','FontName','Times New Roman');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike train
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
axis([xs 0 1])
set(gca,'XTickLabel',[],'XTick',xticks,'YTick',[0 1],'FontName','Times New Roman');
set(gca,'YTick',[0 1],'FontName','Times New Roman');
ylab=ylabel([{'\qquad $n_t$'}; {'($\#$ spikes)'}],'Interpreter','latex','VerticalAlignment','middle','FontName','Times New Roman');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% spike history
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,P.omega*R.h,'Color',gray,'LineWidth',lw)%(R.h-hmin)/(hmax-hmin)
axis([xs round(hmin) round(hmax)])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
ylab=ylabel([{'  $\omega h_t$'}; {'(a.u.)'}],'Interpreter','latex','VerticalAlignment','middle','FontName','Times New Roman');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[round(hmin) round(hmax)],'YTickLabel',{round(hmin); round(hmax)},'FontName','Times New Roman');

% calcium
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
plot(Sim.tvec,R.C,'Color',gray,'LineWidth',lw)%(R.C-cmin)/(cmax-cmin)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
axis([xs round(cmin*100)/100 round(cmax*100)/100])
set(gca,'YTick',[round(cmin*100)/100 round(cmax*100)/100],'YTickLabel',{round(cmin*100)/100; round(cmax*100)/100},'FontName','Times New Roman');
% set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
ylab=ylabel([{'[Ca$^{2+}]_t$'}; {' ($\mu$M)'}],'Interpreter','latex','VerticalAlignment','middle','FontName','Times New Roman');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
% xlab=xlabel('Time (sec)');

% observation state
subplot('Position',[0.166 0.05 0.285 0.08]), cla, hold on, box on
plot(Oind*Sim.dt,O,'.k','LineWidth',1,'markersize',7)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
axis([xs min(O(2:end-6)) max(O(2:end-6))])
xlabel('Fluorescence Observations');
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''},'FontName','Times New Roman');
% title('Observation State','fontsize',tfs)
set(gca,'YTick',[min(O) max(O)],'YTickLabel',{'0'; '1'},'FontName','Times New Roman');
ylab=ylabel([{'\quad $F_t$'}; {'(a.u.)'}],'Interpreter','latex','VerticalAlignment','middle','FontName','Times New Roman');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
xlab=xlabel('Time (sec)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% external stimulus
subplot('Position',[0.62 0.838 0.285 0.08]), cla, hold on
imagesc(Sim.x)
colormap('gray'), 
colorbar([0.92 0.838 0.025 0.08],'YTick',[-0.4  +0.4],'YTickLabel',{'-0.4','+0.4'},'FontName','Times New Roman');
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
axis([xmin xmax 0.5 5.5])
ylab=ylabel({'External      '; 'Stimulus     '});
% ylab=ylabel('$x_t$','Interpreter','latex');
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName','Times New Roman');
% title('Input','fontsize',tfs)

% filtered stimulus
subplot('Position',[0.62 0.73 0.285 0.08]), cla, hold on, box on
plot(Sim.tvec,P.k'*Sim.x,'color','k','linewidth',lw);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color','k','FontName','Times New Roman');
% set(gca,'YAxisLocation','right');
set(gca,'XTick',xticks,'XTickLabel',[],'FontName','Times New Roman');
set(gca,'YTickLabel',[],'FontName','Times New Roman');
axis([xs min(P.k'*Sim.x) max(P.k'*Sim.x)])
ylab=ylabel({'Linearly Filtered'; 'Stimulus    '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName','Times New Roman');

% xlab=xlabel('Time (sec)');

% P(spike)
i=ii+1; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,S.p);
hfill=fill([Sim.tvec Sim.tvec(ind)],([ptiles(1,:) ptiles(2,ind)]-pmin)/(pmax-pmin),ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(R.p-pmin)/(pmax-pmin),'Color',gray,'LineWidth',lw)
plot(Sim.tvec,(M.pbar-pmin)/(pmax-pmin),'linewidth',2,'color',col(2,:))
axis([xs 0 1])
set(gca,'XTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
set(gca,'YTickLabel',[],'FontName','Times New Roman');
% set(gca,'YAxisLocation','right');
% set(gca,'YTick',[0 1])
ylab=ylabel({'Probability   ';'of spiking    '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName','Times New Roman');

% spike train
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
BarVar=M.nbar+M.nvar;
BarVar(BarVar>1)=1;
stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',sw)
stem(Sim.tvec,M.nbar,'Marker','none','Color',col(2,:),'LineWidth',sw)
axis([xs min(R.n) max(R.n)])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
% set(gca,'YAxisLocation','right');
% set(gca,'YTick',[0 1],'YTickLabel',{'0'; '1'})
ylab=ylabel({'Spike Train    '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName','Times New Roman');

% spike history
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,P.omega*S.h);
hfill=fill([Sim.tvec Sim.tvec(ind)],[ptiles(1,:) ptiles(2,ind)],ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,P.omega*R.h,'Color',gray,'LineWidth',lw)%(R.h-hmin)/(hmax-hmin)
plot(Sim.tvec,P.omega*M.hbar,'linewidth',2,'color',col(2,:))
axis([xs round(hmin) round(hmax)])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
% set(gca,'YAxisLocation','right');
% set(gca,'YTick',[round(hmin) round(hmax)],'YTickLabel',{round(hmin); round(hmax)})
ylab=ylabel({'Weighted    ';'Spike History ';'Term        '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName','Times New Roman');

% calcium
i=i+2; subplot(Nrows,Ncols,i), cla, hold on, box on
ptiles = GetPercentiles([.25 .75],S.w_b,S.C);
hfill=fill([Sim.tvec Sim.tvec(ind)],[ptiles(1,:) ptiles(2,ind)],ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,R.C,'Color',gray,'LineWidth',lw)%(R.C-cmin)/(cmax-cmin)
plot(Sim.tvec,M.Cbar,'linewidth',2,'color',col(2,:))
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
axis([xs round(cmin*100)/100 round(cmax*100)/100])
% set(gca,'YAxisLocation','right');
% set(gca,'YTick',[round(cmin*100)/100 round(cmax*100)/100],'YTickLabel',{round(cmin*100)/100; round(cmax*100)/100})
ylab=ylabel({'Calcium     '; 'Concentration '});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName','Times New Roman');
% set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
% xlab=xlabel('Time (sec)');

subplot('Position',[0.62 0.05 0.285 0.08]), cla, hold on, box on
plot(Oind*Sim.dt,O,'.k','LineWidth',1,'markersize',7)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(O(2:end-6)) max(O(2:end-6))])
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''},'FontName','Times New Roman');
% title('Observation State','fontsize',tfs)
ylab=ylabel({'Fluorescence';'Observations'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName','Times New Roman');
axis([xs min(O) max(O)])
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[],'FontName','Times New Roman');
% set(gca,'YAxisLocation','right');
% set(gca,'YTick',[min(O) max(O)],'YTickLabel',{'0'; '1'})
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''},'FontName','Times New Roman');

xlab=xlabel('Time (sec)','FontName','Times New Roman');



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
print -depsc C:\D\Research\liam\SMC_EM_GLM\schem4