ass=t-10:t+10;
ind=find(PHHn==0)

figure(1), clf, hold on
plot(Sim.tvec(ass),S.h(:,ass,:)')
plot(Sim.tvec(ass),S.h(ind,ass,:)','k','linewidth',3)
axis([Sim.tvec(ass(1)) Sim.tvec(ass(end)) min(min(S.h(:,ass,:))) max(max(S.h(:,ass,:)))])

O    = R.O.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);
plot(Sim.tvec(ass),O(ass)','ok','linewidth',2)
plot(Sim.tvec(ass),R.h(ass)','r','linewidth',2)

cmin = min(min(S.C(:,ass)));
cmax = max(max(S.C(:,ass)));

% plot(Sim.tvec(ass),(S.C(:,ass)'-cmin)/(cmax-cmin))
% plot(Sim.tvec(ass),(S.C(ind(2),ass)'-cmin)/(cmax-cmin),'k','linewidth',3)
% plot(Sim.tvec(ass),(R.C(ass)'-cmin)/(cmax-cmin),'r','linewidth',2)
% plot(Sim.tvec(ass),(O(ass)'-cmin)/(cmax-cmin),'ok','linewidth',2)
% stem(repmat(Sim.tvec(ass),numel(ind(2)),1),S.n(ind(2),ass))

% figure(2), clf, hold on
% plot(Sim.tvec(ass),S.C(:,ass)')
% plot(Sim.tvec(ass),S.C(ind,ass)','k','linewidth',3)
% plot(Sim.tvec(ass),R.C(ass)','r','linewidth',2)
% plot(Sim.tvec(ass),O(ass)','ok','linewidth',2)
% stem(repmat(Sim.tvec(ass),numel(ind),1),S.n(ind,ass)*cmax)
% indn=find(S.n(:,t-1));
% % plot(Sim.tvec(ass),S.C(indn,ass)','b','linewidth',3)
% axis([Sim.tvec(ass(1)) Sim.tvec(ass(end)) cmin cmax])
% plot(Sim.tvec(ass),S.C(50,ass)','b','linewidth',3)


