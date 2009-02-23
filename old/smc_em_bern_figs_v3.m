function smc_em_bern_figs_v3(Sim,P,R,S,M)

O=R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);

%%%%%%%% make schematic fig
figure(1)
Nsubs=7;

% for k=1/Sim.dt:Sim.T-1/Sim.dt
%     sumn(k)=sum(R.n(k:k+1/Sim.dt));
% end
% [maxn ind]=max(sumn);
% xs = [ind*Sim.dt-0.75 ind*Sim.dt+0.25];
xmin    = Sim.freq;
xmax    = Sim.T-Sim.freq-1;
xs = [Sim.tvec(xmin) Sim.tvec(xmax)];

subplot(Nsubs,1,1)
if Sim.StimDim==1
    plot(Sim.x,'k')
else
    imagesc(Sim.x), colormap('gray')
    axis([xs/Sim.dt 1 Sim.StimDim])
    % set(gca,'YTick',[])
end
set(gca,'XTick',[])
ylabel('stim')

subplot(Nsubs,1,2)
rateAX  = [xs 0 1];
plot(Sim.tvec,R.p,'k','LineWidth',2)
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
ylabel('rate')
axis(rateAX)

subplot(Nsubs,1,3), hold off
spikemax= max(1,max(R.n(xmin:xmax)));
spikeAX = [xs 0 spikemax];
stem(Sim.tvec,R.n,'Marker','none','Color','k','LineWidth',2)
hold on
plot(Sim.tvec,R.p,'Color',[0.75 0.75 0.75],'LineWidth',1)
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
ylabel('n')
axis(spikeAX)

subplot(Nsubs,1,4), hold off
cmin=min(min(min(R.C(xmin:xmax)),min(M.Cbar(xmin:xmax))),min(O(xmin:xmax)));
cmax=max(max(max(R.C(xmin:xmax)),max(M.Cbar(xmin:xmax))),max(O(xmin:xmax)));
plot(Sim.tvec,(R.C-cmin)/(cmax-cmin),'k','LineWidth',2)
hold on
plot(Sim.tvec,(O-cmin)/(cmax-cmin),'ok','LineWidth',2)
stem(Sim.tvec,R.n,'Marker','none','Color',[0.75 0.75 0.75],'LineWidth',1)
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
ylabel('C')
axis(spikeAX)

subplot(Nsubs,1,5), hold off
plot(Sim.tvec,O-cmin,'.','LineWidth',1,'Color',[0.75 0.75 0.75]);
if Sim.van==true
    Cbar=mean(S.C,1);
else
    Cbar = sum(S.w_f.*S.C,1);
end
Cvar = sum((repmat(Cbar,Sim.N,1)-S.C).^2)/Sim.N;
hold on
plot(Sim.tvec,R.C-cmin,'Color',[0.75 0.75 0.75],'LineWidth',1)
plot(Sim.tvec,Cbar-cmin,'k','LineWidth',2)
plot(Sim.tvec,Cbar+Cvar-cmin,'k','LineWidth',1)
plot(Sim.tvec,Cbar-Cvar-cmin,'k','LineWidth',1)
plot(Sim.tvec,O-cmin,'o','LineWidth',1,'Color',[0.75 0.75 0.75]);
hold off
ylabel('forward')
axis([xs 0 max(Cbar(xmin:xmax)+Cvar(xmin:xmax))-cmin])
% end
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])

subplot(Nsubs,1,6),hold off
plot(Sim.tvec,M.Cbar-cmin,'k','LineWidth',2); hold on
plot(Sim.tvec,R.C-cmin,'Color',[0.75 0.75 0.75],'LineWidth',1)
plot(Sim.tvec,M.Cbar+sqrt(M.Cvar)/2-cmin,'k');
plot(Sim.tvec,M.Cbar-sqrt(M.Cvar)/2-cmin,'k');
% errorbar(Sim.tvec(xmin:xmax),S.Cbar(xmin:xmax),S.Cvar(xmin:xmax),'k','LineWidth',2), hold on,
plot(Sim.tvec,O-cmin,'o','LineWidth',1,'Color',[0.75 0.75 0.75]);
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
ylabel('m_c')
axis([xs 0 max(M.Cbar(xmin:xmax))-cmin])

subplot(Nsubs,1,7), hold off
stem(Sim.tvec,M.nbar,'k','Marker','none','LineWidth',2), hold on
plot(Sim.tvec,M.nbar+sqrt(M.nvar),'k','LineWidth',1);
hold on
stem(Sim.tvec,R.n,'r','Marker','none','LineWidth',1)%'Color',[0.75 0.75 0.75],)
plot(Sim.tvec,(O-cmin)/max(O(xmin:xmax)),'o','LineWidth',1,'Color',[0.75 0.75 0.75]);
% set(gca,'XTickLabel',[])
% set(gca,'YTickLabel',[])
ylabel('m_n')
axis([xs 0 1])

%% %%%%%% plot E[C_t | O], E[C_t | O_{0:k}],C_t, O_t
k=Sim.T;
figure(2), clf, hold on
Sim.tvec=0:Sim.dt:Sim.Nsec-Sim.dt; %time vector
plot(Sim.tvec(1:k),S.C(:,1:k)')
plot(Sim.tvec(1:k),O(1:k),'ok','LineWidth',2)
plot(Sim.tvec(1:k),R.C(1:k),'b','LineWidth',2)
plot(Sim.tvec(1:k),Cbar(1:k),'r','LineWidth',2)
stem(Sim.tvec(1:k),S.n(:,1:k)')
stem(Sim.tvec(1:k),R.n(1:k),'fill')
plot(Sim.tvec(1:k),M.Cbar(1:k),'g','LineWidth',2)
hold off
axis([xs 0 cmax])


%% %%%%%%% plot kernel
% figure, clf, hold on
% if Sim.M>0
%     plot([P.k; P.omega],'-','Color',[0.75 0.75 0.75],'LineWidth',4)
%     for i=1:Sim.iters+1
%         plot([E(i).k; E(i).omega],'k','LineWidth',i*.1)
%     end
%     axis([0 length(P.k)+Sim.M min(min([P.k E(1:i).k]))-1 max(max([P.k E(1:i).k]))+1])
% else
%     plot(P.k,'-','Color',[0.75 0.75 0.75],'LineWidth',4)
%     for i=1:Sim.iters+1
%         plot(E(i).k,'k','LineWidth',i*.1)
%     end
%     axis([0 length(P.k) min(min([P.k E(1:i).k]))-1 max(max([P.k E(1:i).k]))+1])
% end
% xlabel('Parameter','FontSize',24)
% ylabel('Magnitude','FontSize',24)
% legend('Actual Values','Estimated Values','Location','Best')
% title('Linear Filter Estimates','FontSize',24)

% %%%%%%%% plot E[n] vs. n
% figure, clf, hold on
% ntiles=zeros(max(R.n),5);
% clear lens
% for n=0:max(R.n)
%     nums=find(R.n==n);
%     ntiles(n+1,:) = quantile(M.nbar(nums),[.025 .25 .50 .75 .975]);
%     plot(n,ntiles(n+1,2):.01:ntiles(n+1,4),'.k','LineWidth',10)
%     plot(n,ntiles(n+1,1):.01:ntiles(n+1,5),'-k','LineWidth',1)
%     %     plot(n,ntiles(n+1,3),'.b')
% end
% % plot(0:n,0:n,'k')
% xlabel('Estimated n','FontSize',24)
% ylabel('Actual n','FontSize',24)
% title('Box plot','FontSize',24)
