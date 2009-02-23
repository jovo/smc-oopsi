O=R.O.*repmat([NaN*ones(1,Sim.frac-1) 1],1,Sim.K_o);
pf=[2 6];
xs=[100:200];
cmin1=min(S{pf(1)}.C(:,xs));
cmin2=min(S{pf(2)}.C(:,xs));
cmax1=max(S{pf(1)}.C(:,xs));
cmax2=max(S{pf(2)}.C(:,xs));
%%
k=Sim.K;
figure, clf, 
subplot(2,1,1), hold on
plot(Sim.tvec(xs),(S{pf(1)}.C(:,xs)'-cmin1)/(cmax1-cmin1),'Color',[0.75 0.75 0.75])
plot(Sim.tvec(xs),(O(xs)-cmin1)/(cmax1-cmin1),'ok','LineWidth',2)
plot(Sim.tvec(xs),(R.C(xs)-cmin1)/(cmax1-cmin1),'k','LineWidth',2)

subplot(2,1,2), hold on
plot(Sim.tvec(xs),(S{pf(2)}.C(:,xs)'-cmin2)/(cmax2-cmin2),'Color',[0.75 0.75 0.75])
plot(Sim.tvec(xs),(O(xs)-cmin2)/(cmax2-cmin2),'ok','LineWidth',2)
plot(Sim.tvec(xs),(R.C(xs)-cmin2)/(cmax2-cmin2),'k','LineWidth',2)


%%
% figure
% Sim.tvec=0:Sim.dt:Sim.Nsec-Sim.dt; %time vector
% Nsubs=6;
% 
% for k=1:Sim.K-1/Sim.dt
%     sumn(k)=sum(R.I(k:k+1/Sim.dt));
% end
% [maxn ind]=max(sumn);
% xs = [ind*Sim.dt ind*Sim.dt+1];
% 
% subplot(Nsubs,1,1)
% if Sim.KernelSize==1
%     plot(Sim.x,'k')
% else
%     imagesc(Sim.x), colormap('gray')
%     axis([xs/Sim.dt 1 Sim.KernelSize])
%     % set(gca,'YTick',[])
% end
% set(gca,'XTick',[])
% ylabel('stim')
% 
% subplot(Nsubs,1,2)
% xmin    = find(Sim.tvec>xs(1),1);
% xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
% ratemax = max(R.p(xmin:xmax));
% rateAX  = [xs 0 ratemax];
% plot(Sim.tvec,R.p,'k','LineWidth',2)
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% ylabel('rate')
% axis(rateAX)
% 
% subplot(Nsubs,1,3), hold off
% spikeAX = [xs 0 1];
% stem(Sim.tvec,R.I,'Marker','none','Color','k','LineWidth',2)
% hold on
% plot(Sim.tvec,R.p*Sim.dt/(ratemax),'Color',[0.75 0.75 0.75],'LineWidth',1)
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% ylabel('n')
% axis(spikeAX)
% 
% subplot(Nsubs,1,4), hold off
% cmin=min(min(min(R.C(xmin:xmax)),min(M(1,1).Cbar(xmin:xmax))),min(O(xmin:xmax)));
% cmax=max(max(max(R.C(xmin:xmax)),max(M(1,1).Cbar(xmin:xmax))),max(O(xmin:xmax)));
% plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'k','LineWidth',2)
% hold on
% stem(Sim.tvec,R.I,'Marker','none','Color',[0.75 0.75 0.75],'LineWidth',1)
% plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',2)
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% ylabel('C')
% axis(spikeAX)
% 
% subplot(Nsubs,1,5),hold off, hold on
% plot(Sim.tvec,M(1,1).Cbar-cmin,'r','LineWidth',2); 
% plot(Sim.tvec,M(1,1).Cbar+sqrt(M(1,1).Cvar)/2-cmin,'r');
% plot(Sim.tvec,M(1,1).Cbar-sqrt(M(1,1).Cvar)/2-cmin,'r');
% 
% plot(Sim.tvec,M(3,1).Cbar-cmin,'b','LineWidth',2); hold on
% plot(Sim.tvec,M(3,1).Cbar+sqrt(M(3,1).Cvar)/2-cmin,'b');
% plot(Sim.tvec,M(3,1).Cbar-sqrt(M(3,1).Cvar)/2-cmin,'b');
% 
% plot(Sim.tvec,(R.C-cmin),'k','LineWidth',2)
% plot(Sim.tvec,O-cmin,'.','LineWidth',1,'Color',[0.75 0.75 0.75]);
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% ylabel('m_c')
% axis([xs 0 max(M(1,1).Cbar(xmin:xmax))-cmin])
% 
% subplot(Nsubs,1,6), hold off, hold on
% stem(Sim.tvec,M(1,1).Ibar,'r','Marker','none','LineWidth',2)
% plot(Sim.tvec,M(1,1).Ibar+sqrt(M(1,1).Ivar),'r','LineWidth',1);
% 
% stem(Sim.tvec,M(3,1).Ibar,'b','Marker','none','LineWidth',2)
% plot(Sim.tvec,M(3,1).Ibar+sqrt(M(3,1).Ivar),'b','LineWidth',1);
% 
% stem(Sim.tvec,R.I,'k','LineWidth',1)
% plot(Sim.tvec,(O-cmin)/max(O(xmin:xmax)),'.','LineWidth',1,'Color',[0.75 0.75 0.75]);
% % set(gca,'XTickLabel',[])
% % set(gca,'YTickLabel',[])
% ylabel('m_n')
% axis([xs 0 1])

%% 
figure, hold on
plot(Sim.tvec,R.C,'k','LineWidth',2)
plot(Sim.tvec,O,'ok')
plot(Sim.tvec,M(1,1).Cbar,'r')
plot(Sim.tvec,M(5,1).Cbar,'g')

% figure, imagesc(S.w_b)