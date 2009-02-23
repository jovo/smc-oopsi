function GetFig_Array2(Sim,R,M,Os,freq)

%% preset stuff for fig
figure(3), clf, cla                 %clear the fig
set(gcf, 'color', 'w');             %make background color white

Nrows   = length(M);                %rows are different observation noise
Ncols   = length(freq);             %columns are different sampling frequencies

gray    = [0.75 0.75 0.75];         %define gray
col     = [0 0 1; 0 .5 0; 1 0 0];   %define colors for mean
ccol    = col+.8; ccol(ccol>1)=1;   %define colors for std
ind     = Sim.T:-1:1;               %inverse indices for 'fill' function

tl   = [.05 0.25];                  %tick length
yfs  = 16;                          %ylabel font size
xfs  = 14;                          %xlabel font size

xs      = [Sim.tvec(Sim.freq) Sim.tvec(end-Sim.freq*2)];%the limits of x-axis on which to plot things
xmin    = find(Sim.tvec>xs(1),1);   %find index of first time step
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;%find index of last time step
xind    = xmin:xmax;                %indices of x-axis

for i=1:Ncols                       %find effective min and max for each simulation
    for n=1:Nrows
        cmin(i,n)=min(M(n,i).Cbar(xind)-sqrt(M(n,i).Cvar(xind)));
        cmax(i,n)=max(M(n,i).Cbar(xind)+sqrt(M(n,i).Cvar(xind)));
    end
end
cmin=min(cmin(:));                  %use cmin and cmax to normalize...
cmax=max(cmax(:));                  %all the calcium plots to go between 0 and 1


%% make fig
for i=1:Ncols
    Sim.freq    = freq(i);                                                              %set frequency
    Sim.T_o     = Sim.T/Sim.freq;                                                       %number of observations
    for n=1:Nrows
        O  = Os{n,i}.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.T_o);                     %recall observations
        j=(Nrows-n)*(Ncols)+i;                                                          %index of subplot
        subplot(Nrows,Ncols,j), cla, hold on

        % plot calcium stuff
        hfill=fill([Sim.tvec Sim.tvec(ind)],...                                         %fill space between +/- std of mean
            ([M(n,i).Cbar-sqrt(M(n,i).Cvar) M(n,i).Cbar(ind)+sqrt(M(n,i).Cvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
        set(hfill,'edgecolor',ccol(2,:))                                                %make edge of fill color same as fill color
        plot(Sim.tvec,(M(n,i).Cbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))  %plot mean
        plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',1,'markersize',4)         %plot observations
        plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)              %plot true calcium
        plot(Sim.tvec,ones(size(Sim.tvec)),'k')                                         %plot a line dividing calcium and spikes

        % plot spike stuff
        stem(Sim.tvec,R.n,'Marker','none','Color',gray,'LineWidth',1)                   %plot actual spikes
        BarVar=M(n,i).nbar+M(n,i).nvar;                                                 %make sure var of spikes doesn't get larger than 1
        BarVar(BarVar>1)=1;
        stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',2)           %plot variance of spikes
        stem(Sim.tvec,M(n,i).nbar,'Marker','none','Color',col(2,:),'LineWidth',2)       %plot mean of spikes
        axis([xs 0 2])

        % set labels and stuch
        if n==1, set(get(gca,'XLabel'),'String',[num2str(freq(i)/Sim.dt) ' Hz'],'fontsize',xfs), end
%         if n==1
%             if i==1, set(get(gca,'XLabel'),'String','200 Hz','fontsize',xfs);
%             elseif i==2, set(get(gca,'XLabel'),'String','40 Hz','fontsize',xfs);
%             elseif i==3, set(get(gca,'XLabel'),'String','20 Hz','fontsize',xfs);
%             elseif i==4, set(get(gca,'XLabel'),'String','10 Hz','fontsize',xfs);
%             elseif i==5, set(get(gca,'XLabel'),'String','5 Hz','fontsize',xfs);
%             end
%         end

        %         if j==1+(Nrows-1)*Ncols, set(get(gca,'YLabel'),'String',texlabel([num2str(freq(i)) ' sigma_o']),'fontsize',yfs);   end
        %         if j==1+(Nrows-1)*Ncols, set(get(gca,'YLabel'),'String',texlabel([num2str(freq(i)) ' sigma_o']),'fontsize',yfs);   end
        %         if j==1+(Nrows-2)*Ncols, set(get(gca,'YLabel'),'String',texlabel([num2str(freq(i)) ' sigma_o']),'fontsize',yfs);  end
        %         if j==1+(Nrows-3)*Ncols, set(get(gca,'YLabel'),'String',texlabel('5 sigma_o'),'fontsize',yfs);  end
        %         if j==1+(Nrows-4)*Ncols, set(get(gca,'YLabel'),'String',texlabel('10 sigma_o'),'fontsize',yfs);  end
        %         if j==1+(Nrows-5)*Ncols, set(get(gca,'YLabel'),'String',texlabel('20 sigma_o'),'fontsize',yfs);  end

        set(gca,'TickLength',tl,'XTick',xs(1):.1:xs(2))
        set(gca,'XTickLabel',[],'YTickLabel',[]),
    end
end

%% add some stuff to make fig pretty
annotation('arrow',[0.07 0.07],[0.4084 0.628],'linewidth',2);   	%yaxis arrow  %[x(1) x(2)], [y(1) y(2)]
annotation('arrow',[0.6219 0.4156],[0.07 0.07],'linewidth',1.5);    %xaxis arrow

text('Interpreter','tex',...
    'String','Increasing Sampling Frequency',...
    'Position',[-2.355 -12.33 17.32],...
    'FontSize',14,...
    'FontName','Helvetica')

text('Interpreter','tex',...
    'String','Increasing Observation Noise',...
    'Position',[-4.588 -7.46 17.32],...
    'FontSize',14,...
    'FontName','Helvetica',...
    'Rotation',90)


%% print fig to epsc
fig=figure(3);
bgr=1*[8.5 8.5];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\array