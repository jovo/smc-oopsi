fig=figure(101); clf,
gray  = [.75 .75 .75];              % define gray color
col   = [1 0 0; 0.2 0.2 1];            % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter = 'none';                     % interpreter for axis labels
xlims = [4 Sim.T-2];                % xmin and xmax for current plot
fs=12;                              % font size
fn=['Helvetica'];                   % font name
in='tex';                           % interpreter
yfs=12;                             % ylabel font size
xfs=12;                             % xlabel font size
ms=20;                              % marker size for real spike
ms2=3.5;                              % marker size for real spike
sw=3.5;                             % spike width
lw=2;                               % line width
tl=[.05 0.25];                      %tick length
i=0;
xtick=[15 30 45];%vec(end):vec(end):Sim.T;
yshift = .05;

for sub=1:length(vec)
    for nois=1:length(vec)
        subplot('Position',[(nois-1)*0.18+0.15 (length(vec)-sub)*0.20+.12 0.17 0.18]) %[left bottom width height]
        hold on     
        obs=vec(sub):vec(sub):Sim.T;
        plot(obs,z1(I{sub,nois}.F(obs))+1+yshift,'.k','LineWidth',lw,'MarkerSize',10);
        %         stem(spt,Sim.n(spt),'Marker','none','MarkerSize',ms,'LineWidth',sw,'Color',gray)
        BarVar=I{sub,nois}.M.nbar+I{sub,nois}.M.nvar; BarVar(BarVar>1)=1;
        spts=find(BarVar>1e-3);
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',gray);
%         bar(spts,BarVar(spts),'EdgeColor',gray,'FaceColor','w');%,'Marker','none','LineWidth',sw,'Color',gray);        
        spts=find(I{sub,nois}.M.nbar>1e-3);
        stem(spts,I{sub,nois}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color','k')
        stem(spt,1.1*Sim.n(spt),'Marker','v','MarkerSize',ms2,'LineStyle','none','MarkerFaceColor',gray,'MarkerEdgeColor',gray)
        axis([xlims 0 2+yshift])
        axis('tight')
        box('off')
        set(gca,'YTick',0:2,'YTickLabel',[],'XTick',xtick,'XTickLabel',[],'FontName',fn,'TickLength',tl)
        
        plot(ones(1,Sim.T),'k')
        
        if nois==1,
            xlab=[num2str(vec(sub)*Sim.dt*1000) ' ms'];
            set(get(gca,'YLabel'),'String',texlabel(xlab),'fontsize',yfs,'color','k','interpreter',in...
                ,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','FontName',fn)
        end
        if sub==4,
            ylab=[num2str(vec(nois)) ' \sigma_F'];
            set(get(gca,'XLabel'),'String',ylab,'fontsize',xfs,'interpreter','tex')%,'FontName',fn),
        end
        if sub==4 && nois==1
            set(gca,'YTick',0:2,'YTickLabel',[],'XTick',xtick,'XTickLabel',Sim.tvec(xtick),'FontName',fn)
            set(get(gca,'XLabel'),'String','Time (sec)','FontSize',fs,'Interpreter',in,'FontName',fn); %xlabel('Time (sec)','FontSize',fs)
        end
    end
end

xpos=0.05;
ratepos=0.05;
annotation('arrow',[0.4 0.59],[xpos xpos],'LineWidth',lw,'HeadStyle','plain');
% annotation('arrow',[0 0],[0 0],'Color',[1 1 1]);
annotation('arrow',[ratepos ratepos],[0.39 0.62],'LineWidth',lw,'HeadStyle','plain');

text('Position',[-170 6.3],...
    'Interpreter',in,...
    'String','Increasing Frame Rate',...
    'FontSize',fs,...
    'FontName',fn,...
    'Rotation',90,...
    'HorizontalAlignment','right',...
    'verticalalignment','middle')

text('Position',[-85 -1.12],...
    'Interpreter',in,...
    'String','Increasing Observation Noise',...
    'FontSize',fs,...
    'FontName',fn,...
    'Rotation',0)

% print fig
wh=[7 5];   %width and height
set(gcf,'PaperSize',wh);
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-depsc','ArraySim_bw')