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
%     subplot(3,2,siz(1)-nn+1), hold on
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
%     
%     if nn==5
%         title('Linear Filter','fontsize',tfs)
%     end
%     
%     subplot(3,2,siz(1)-nn+2), hold on
%     
%     if nn==5
%         title('Calcium Parameters','fontsize',tfs)
%     end
%     
% end
% fig=figure(8);
% bgr=[6 6];
% set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
% print -depsc C:\D\Research\liam\SMC_EM_GLM\mstep1



%%
figure(9), clf, hold on
load Ca_params2
nSpikes     = 50.*2.^[0:4];
siz=size(E);
params=[P.tau_c; P.A; P.C_0; P.sigma_c];
for nn=1:3%:siz(1)
    clear ca_est 
    if nn==3
        for kk=1:4
            ca_est(kk,1)=E(nn,kk).tau_c;
            ca_est(kk,2)=E(nn,kk).A;
            ca_est(kk,3)=E(nn,kk).C_0;
            ca_est(kk,4)=E(nn,kk).sigma_c;
        end
    else
        for kk=1:siz(2)
            ca_est(kk,1)=E(nn,kk).tau_c;
            ca_est(kk,2)=E(nn,kk).A;
            ca_est(kk,3)=E(nn,kk).C_0;
            ca_est(kk,4)=E(nn,kk).sigma_c;
        end
    end
    meanca_est=mean(ca_est);
    stdca_est=sqrt(var(ca_est));
    subplot(3,2,siz(1)-nn+1), cla, hold on
    %     subplot('Position', [.1 b+h*(nn-1) .43 h]), hold on
    
    if nn==5
        title('Linear Filter','fontsize',tfs)
    end
    
    subplot(3,2,siz(1)-nn+2), hold on
    plot(params,'o','linewidth',2,'color',gray)
    errorbar(meanca_est,stdca_est,'.','color','k','linewidth',2)
    set(gca,'XTickLabel',[],'YTickLabel',[])
    axis([.9 4.1 min(0,min(meanca_est-stdca_est)) max(max(params),max(meanca_est-stdca_est))])
    ylab=ylabel({nSpikes(nn);'Spikes'});
    set(ylab,'Rotation',0,...
        'HorizontalAlignment','right',...
        'verticalalignment','middle',...
        'fontsize',yfs)
    
    if nn==5
        title('Calcium Parameters','fontsize',tfs)
    end
    
end
fig=figure(8);
bgr=[6 6];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\mstep1


% %%
% ylab=ylabel('Magnitude (a.u.)');
% set(ylab,'fontsize',yfs)
% title('Linear Filter','fontsize',tfs)
% 
% subplot('Position', [.60 bl .35 .8]), cla, hold on
% thetastar=[P.tau_c; P.A; P.C_0]; hold on
% % plot(thetastar,'color',gray,'linewidth',2)
% plot(thetastar./thetastar,'o','color',gray,'linewidth',2)
% for kk=2:length(E)
%     theta(kk,:)=[E(kk).tau_c; E(kk).A; E(kk).C_0];
% end
% meantheta=mean(theta);
% stdtheta=sqrt(var(theta));
% errorbar(meantheta./thetastar',stdtheta,'color','k','linewidth',2)
% axis([.9 3.1 min(0,min(meantheta./thetastar'-stdtheta)) max(meantheta./thetastar'+stdtheta)])
% 
% tit='Calcium Parameters';
% title(tit,'fontsize',tfs)
% 
% ypos = -0.25;
% text('Position',[0.9517 ypos 17.32],...
%     'Interpreter','tex',...
%     'String','{\tau}_{c}',...
%     'FontSize',yfs);
% text('Position',[1.98 ypos 17.32],...
%     'Interpreter','tex',...
%     'String','A',...
%     'FontSize',yfs);
% text('Position',[2.8 ypos 17.32],...
%     'Interpreter','tex',...
%     'String','[Ca^{2+}]_0',...
%     'FontSize',yfs);
% set(gca,'XTick',[1 2 3],'XTickLabel',[],'YTick',[1 2 3])
% 
% ylab=ylabel('Error (%/100)');
% set(ylab,'fontsize',yfs)
% 
% %
% fig=figure(8);
% bgr=[6 2];
% set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
% print -depsc C:\D\Research\liam\SMC_EM_GLM\mstep

% %%
% ylab=ylabel('Magnitude (a.u.)');
% set(ylab,'fontsize',yfs)
% title('Linear Filter','fontsize',tfs)
% 
% subplot('Position', [.60 bl .35 .8]), cla, hold on
% thetastar=[P.tau_c; P.A; P.C_0]; hold on
% % plot(thetastar,'color',gray,'linewidth',2)
% plot(thetastar./thetastar,'o','color',gray,'linewidth',2)
% for kk=2:length(E)
%     theta(kk,:)=[E(kk).tau_c; E(kk).A; E(kk).C_0];
% end
% meantheta=mean(theta);
% stdtheta=sqrt(var(theta));
% errorbar(meantheta./thetastar',stdtheta,'color','k','linewidth',2)
% axis([.9 3.1 min(0,min(meantheta./thetastar'-stdtheta)) max(meantheta./thetastar'+stdtheta)])
% 
% tit='Calcium Parameters';
% title(tit,'fontsize',tfs)
% 
% ypos = -0.25;
% text('Position',[0.9517 ypos 17.32],...
%     'Interpreter','tex',...
%     'String','{\tau}_{c}',...
%     'FontSize',yfs);
% text('Position',[1.98 ypos 17.32],...
%     'Interpreter','tex',...
%     'String','A',...
%     'FontSize',yfs);
% text('Position',[2.8 ypos 17.32],...
%     'Interpreter','tex',...
%     'String','[Ca^{2+}]_0',...
%     'FontSize',yfs);
% set(gca,'XTick',[1 2 3],'XTickLabel',[],'YTick',[1 2 3])
% 
% ylab=ylabel('Error (%/100)');
% set(ylab,'fontsize',yfs)
% 
% %
% fig=figure(8);
% bgr=[6 2];
% set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
% print -depsc C:\D\Research\liam\SMC_EM_GLM\mstep