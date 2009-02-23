function GetFig_Array1(Sim,R,M,Os,freq)

Ncols=length(freq);
Nrows=length(M);
xs      = [Sim.tvec(Sim.freq) Sim.tvec(end-Sim.freq*2)];
xmin    = find(Sim.tvec>xs(1),1);
xmax    = find(Sim.tvec>xs(2),1)+Sim.dt;
xind    = [xmin:xmax];

gray=[0.75 0.75 0.75];

col=[0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];
ccol=col+.8; ccol(ccol>1)=1;
ind=Sim.K:-1:1;

for i=1:Ncols
    for n=1:Nrows
        cmin(i,n)=min(M(n,i).bCbar-sqrt(M(n,i).bCvar));
        cmax(i,n)=max(M(n,i).bCbar+sqrt(M(n,i).bCvar));
    end
end
cmin=min(cmin(:));
cmax=max(cmax(:));

tl   = [.05 0.25];
yfs  = 16;
xfs  = 14;

figure(3), cla, clf,
for i=1:Ncols
    Sim.freq    = freq(i);
    Sim.K_o     = Sim.K/Sim.freq;                %number of observations
    for n=1:Nrows
        O  = Os{n,i}.*repmat([NaN*ones(1,Sim.freq-1) 1],1,Sim.K_o);
        j=(Nrows-n)*(Ncols)+i;
        subplot(Ncols,Nrows,j), cla, hold on

        h=fill([Sim.tvec Sim.tvec(ind)],([M(n,i).bCbar-sqrt(M(n,i).bCvar) M(n,i).bCbar(ind)+sqrt(M(n,i).bCvar(ind))]-cmin)/(cmax-cmin)+1,ccol(2,:));
        set(h,'edgecolor',ccol(2,:))
        plot(Sim.tvec,(M(n,i).bCbar-cmin)/(cmax-cmin)+1,'linewidth',2,'color',col(2,:))
        plot(Sim.tvec,(O-cmin)/(cmax-cmin)+1,'ok','LineWidth',1,'markersize',4)
        plot(Sim.tvec,(R.C-cmin)/(cmax-cmin)+1,'color',gray,'LineWidth',1)
        plot(Sim.tvec,ones(size(Sim.tvec)),'k')

        stem(Sim.tvec,R.I,'Marker','none','Color',gray,'LineWidth',1)
        BarVar=M(n,i).bIbar+M(n,i).bIvar;
        BarVar(BarVar>1)=1;
        stem(Sim.tvec,BarVar,'Marker','none','Color',ccol(2,:),'LineWidth',2)
        stem(Sim.tvec,M(n,i).bIbar,'Marker','none','Color',col(2,:),'LineWidth',2)
        axis([xs 0 2])

        if j==Ncols*(Nrows-1)+(Ncols-4), set(get(gca,'XLabel'),'String','200 Hz','fontsize',xfs);    end
        if j==Ncols*(Nrows-1)+(Ncols-3), set(get(gca,'XLabel'),'String','40 Hz','fontsize',xfs);   end
        if j==Ncols*(Nrows-1)+(Ncols-2), set(get(gca,'XLabel'),'String','20 Hz','fontsize',xfs);   end
        if j==Ncols*(Nrows-1)+(Ncols-1), set(get(gca,'XLabel'),'String','10 Hz','fontsize',xfs);  end
        if j==Ncols*(Nrows-1)+(Ncols-0), set(get(gca,'XLabel'),'String','5 Hz','fontsize',xfs);  end
        
        if j==1+(Nrows-1)*Ncols, set(get(gca,'YLabel'),'String',texlabel('0'),'fontsize',yfs);   end
        if j==1+(Nrows-2)*Ncols, set(get(gca,'YLabel'),'String',texlabel('sigma_o'),'fontsize',yfs);  end
        if j==1+(Nrows-3)*Ncols, set(get(gca,'YLabel'),'String',texlabel('2 sigma_o'),'fontsize',yfs);  end
        if j==1+(Nrows-4)*Ncols, set(get(gca,'YLabel'),'String',texlabel('5 sigma_o'),'fontsize',yfs);  end
        if j==1+(Nrows-5)*Ncols, set(get(gca,'YLabel'),'String',texlabel('10 sigma_o'),'fontsize',yfs);  end

        set(gca,'TickLength',tl,'XTick',[xs(1):.1:xs(2)])
        set(gca,'XTickLabel',[])%[xs(1):.1:xs(2)]-xs(1))

        set(gca,'YTickLabel',[]),
    end
end

annotation('arrow',[0.07 0.07],[0.4084 0.628],'linewidth',2);   	%yaxis arros  %[x(1) x(2)], [y(1) y(2)]
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

        
fig=figure(3);
bgr=1*[8.5 8.5];
set(fig,'PaperPosition',[0 11-bgr(2) bgr]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\array