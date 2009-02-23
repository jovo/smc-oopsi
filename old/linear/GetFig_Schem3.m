function GetFig_Schem3(Sim,R,S,M)

%% initialize stuff
%let O be only observations at sample times
O           = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
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
cmin    = min(min(min(R.C(xind),min(M.Cbar(xind)-sqrt(M.Cvar(xind))))),min(O(2:end-1)));
cmax    = max(max(max(R.C(xind),max(M.Cbar(xind)+sqrt(M.Cvar(xind))))),max(O(2:end-1)));

gray    = [0.75 0.75 0.75];
col     = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol    = col+.8; ccol(ccol>1)=1;
ind     = Sim.T:-1:1;
sw      = 4.5;
fs      = 14;

figure(1), clf, Nsubs=4;
set(gcf, 'color', 'w');

%% external stimulus
i=1; subplot(Nsubs,1,i), cla, hold on
plot(Sim.tvec,Sim.x,'k','LineWidth',2)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis([xs min(Sim.x(xind)) max(Sim.x(xind))+1.])
ylab=ylabel({'External';'Stimulus'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% hidden states
i=i+1; subplot(Nsubs,1,i), cla, hold on
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'Color',gray,'LineWidth',2)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis(hidAX)
ylab=ylabel({'True [Ca^2^+]';'and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',gray)

%% observation state
i=i+1; subplot(Nsubs,1,i), cla, hold on
plot(Oind*Sim.dt,(O-cmin)/(cmax-cmin),'.k','LineWidth',1,'markersize',7)
set(gca,'XTickLabel',[],'YTickLabel',[],'XTick',xticks,'YTick',[])
axis(spikeAX)
ylab=ylabel({'Observations'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

%% inference
i=i+1; subplot(Nsubs,1,i), cla, hold on,

% calcium
hfill=fill([Sim.tvec Sim.tvec(ind)],([M.Cbar+sqrt(M.Cvar) M.Cbar(ind)-sqrt(M.Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
set(hfill,'edgecolor',ccol(2,:))
plot(Sim.tvec,(M.Cbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
plot(Sim.tvec,ones(size(Sim.tvec)),'k')

% spikes
stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',sw)
BarVar=M.nbar+M.nvar;
BarVar(BarVar>1)=1;
stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',sw)
stem(Sim.tvec,M.nbar,'Marker','none','Color',col(2,:),'LineWidth',sw)

% set stuff
set(gca,'YTick',[]), set(gca,'YTickLabel',[],'XTick',xticks),
set(gca,'XTickLabel',{xticks(1);'';xticks(3);'';xticks(5);'';xticks(7);'';xticks(9);'';xticks(11);'';xticks(13);''})
axis(hidAX)
ylab=ylabel({'Inferred [Ca^2^+]';'and';'Spike Train'});
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','color',col(2,:))
xlab=xlabel('Time (sec)');

% label subplots
% text('String','(A)',...
%     'Position',[0.03 10.71 17.32]);%,...'FontSize',fs);
% 
% text('String','(B)',...
%     'Position',[0.03 8.015 17.32]);%,...'FontSize',fs);
% 
% text('String','(C)',...
%     'Position',[0.03 4.985 17.32]);%,...'FontSize',fs);
% 
% text('String','(D)',...
%     'Position',[0.03 2.348 17.32]);%,...'FontSize',fs);

% figure(12), clf, hold on
% plot(S.C'+cmax), 
% plot(R.C+cmax,'k','linewidth',2), 
% plot(M.Cbar+cmax,'r','linewidth',2),     
% plot(Oind,O+cmax,'.k','LineWidth',1,'markersize',7)      %plot observations



%% make eps
fig=figure(1);
wh=[7 3.5];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\bernoulli\schem