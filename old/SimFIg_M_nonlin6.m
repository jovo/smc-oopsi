%%
clear, clc
gray    = [0.75 0.75 0.75];
yfs     = 16;
tfs     = 18;

%%
% load Rate_params
% figure(8), clf, hold on
% b       = .05;
% h       = .3;
% siz=size(E);
% for nn=1:2:siz(1)
%     clear kvec theta
%     if nn==5
%         for kk=1:2
%             kvec(kk,:)=E(nn,kk).k;
%         end
%     else
%         for kk=1:siz(2)
%             kvec(kk,:)=E(nn,kk).k;
%         end
%     end
%     meankvec=mean(kvec);
%     stdkvec=sqrt(var(kvec));
%     if nn==1, subplot(2,2,3), cla, hold on
%     elseif nn==3, subplot(2,2,1), cla, hold on,
%     elseif nn==5, subplot(2,2,2), cla, hold on
%     end
%     %     subplot('Position', [.1 b+h*(nn-1) .43 h]), hold on
%     plot(P.k,'linewidth',2,'color',gray)
%     errorbar(meankvec,stdkvec,'color','k','linewidth',2)
%     set(gca,'XTickLabel',[],'YTickLabel',[])
%     axis([.9 Sim.StimDim+.1 min(min(P.k),min(meankvec-stdkvec)) max(max(P.k),max(meankvec-stdkvec))])
%     ylab=ylabel({nSpikes(nn);'Spikes'});
%     set(ylab,'Rotation',0,...
%         'HorizontalAlignment','right',...
%         'verticalalignment','middle',...
%         'fontsize',yfs)    
% end

load Ca_params2
params=[P.tau_c; P.A; P.C_0; P.sigma_c];
subplot(2,2,4), cla, hold on
% title('Calcium Parameters','fontsize',tfs)
nn=1;        
for kk=1:siz(2)
    ca_est(kk,1)=E(nn,kk).tau_c;
    ca_est(kk,2)=E(nn,kk).A;
    ca_est(kk,3)=E(nn,kk).C_0;
    ca_est(kk,4)=E(nn,kk).sigma_c;
end
meanca_est=mean(ca_est);
stdca_est=sqrt(var(ca_est));
plot(params,'o','linewidth',2,'color',gray)
errorbar(meanca_est,stdca_est,'+','markersize',7,'color','k','linewidth',2)
% [hx,hy] = format_ticks(gca,{'$\tau_c$','$A$','[Ca$^{2+}]_0$','$\sigma_c$'},[],[1,2,3,4],'FontSize',16);
[hx,hy] = format_ticks(gca,...
    {'$\tau_c$','$A$','[Ca$^{2+}]_0$','$\sigma_c$'},...
    '',[1,2,3,4],[],0,45,[],...
    'FontSize',16,'FontWeight','Bold');
set(gca,'YTick',[.25 .5],'YTickLabel',[])
axis([.9 3.1 min(0,min(meanca_est-stdca_est)) max(max(params),max(meanca_est-stdca_est))])
    
fig=figure(8);
bgr=[6 6];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\mstep2