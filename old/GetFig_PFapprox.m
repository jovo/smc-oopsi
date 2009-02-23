function GetFig_PFapprox(O)

gray= [0.75 0.75 0.75];
col = [0 0 1; 0 .5 0; 1 0 0; 0 1 1; 1 0 1; 1 .5 0; 1 .5 1];

t   = -10:.01:10;
mu1 = -5;
var1= 4;
mu2 = mu1+2.2*var1;
var2= 10;
n   = 50;
rs  = .011;
ms  = 28;
bs  = -.01;
fs  = 18;
xfs = fs;
yfs = fs;
nfs = fs-4;

gaussmix = 1/sqrt(2*pi*var1)*exp(-(mu1-t).^2/var1);
gaussmix = gaussmix+1/sqrt(2*pi*var2)*exp(-(mu2-t).^2/var2);

row{1} = -8:8;
row{2} = row{1}; row{2}([1 8])=[];
row{3} = row{2}; row{3}([6 end])=[];
row{4} = row{3}; row{4}([6])=[];
row{5} = row{4}; row{5}([ end])=[];
row{6} = row{5}; row{6}([6])=[];
row{7} = row{6}; row{7}([1])=[];
row{8} = row{7}; row{8}([4 end])=[];
row{9} = row{8}; row{9}([4])=[];
row{10} = row{9}; row{10}([])=[];
row{11} = row{10}; row{11}([end])=[];
row{12} = row{11}; row{12}([end-1])=[];
row{13} = row{12}; row{13}([end])=[];
row{14} = row{13};
row{15} = row{14}; row{15}([1])=[];
row{16} = row{15}; row{16}([end])=[];
row{17} = row{16};
row{18} = row{17};

figure(5), clf, subplot(1,2,1), hold on,
plot(-gaussmix,t,'k','linewidth',2)

for i = 1:length(row)
    plot(-rs*i*ones(size(row{i}))-bs,row{i},'.','markersize',ms,'color',gray)
end
axis([-max(gaussmix)*1.01 0 min(t) max(t)])

set(gca,'Xtick',[-.2:.1:0],'XTickLabel',[.2:-.1:0],'fontsize',nfs)
set(gca,'Ytick',[mu1 mu2],'YTickLabel',[1 2],'fontsize',nfs)
xlab=xlabel({'Probability'});
set(xlab,'fontsize',xfs)

ylab=ylabel({'Calcium (a.u.)'});
set(ylab,'fontsize',yfs)

subplot(1,2,2), hold on
x       = 15:-.01:-15;
len     = length(O.sig2_o);
O.mu_o(O.mu_o==0)=NaN;

gaussA   = zeros(len,length(x));
gaussB  = zeros(len,length(x));
for u=len:-1:1
    for m=1:2^(len-u)
        gaussA(u,:) = gaussA(u,:) + O.p_o(m,u)*1/sqrt(2*pi*O.sig2_o(u)) * exp(-(x-O.mu_o(m,u)).^2/O.sig2_o(u));
    end
    for m=1:len
        gaussB(u,:) = gaussB(u,:) + O.p(m,u)*1/sqrt(2*pi*O.sig2_o(u)) * exp(-(x-O.mu(m,u)).^2/O.sig2_o(u));
    end
    gaussA(u,:) = gaussA(u,:)/sum(gaussA(u,:));
    gaussB(u,:) = gaussB(u,:)/sum(gaussB(u,:));
end
mxA=max(gaussA(:));
mxB=max(gaussB(:));

for u=len:-1:1
    plot(-gaussA(u,:)/mxB+u,x+len,'k','linewidth',2)
    plot(-gaussB(u,:)/mxB+u,x+len,'color',gray,'linewidth',1)
       for m=1:2^(len-u)
            [foo ind]=find(x<O.mu_o(m,u),1);
            plot(u-gaussA(u,ind)/mxA,len+O.mu_o(m,u),'+k','markersize',10,'linewidth',1)
       end
       for m=1:len-u+1
            [foo ind]=find(x<O.mu(m,u),1);
            plot(u-gaussB(u,ind)/mxB,len+O.mu(m,u),'xk','markersize',10,'linewidth',1,'color',gray)
       end
end
maxmean=max(O.mu_o(:));
minmean=min(O.mu_o(:));
maxvar=4*sqrt(max(O.sig2_o(:)));
axis([0 len minmean-maxvar+len maxmean+maxvar+len])

set(gca,'XTick',[0:5],'XTickLabel',{'','v-4','v-3','v-2','v-1','v'},'fontsize',nfs)
set(gca,'YTick',[1:len]+O.mu_o(1,len),'YTickLabel',[1:len],'fontsize',nfs)
xlab=xlabel({'Time Step'});
set(xlab,'fontsize',xfs)


%% make eps
fig=figure(5);
wh=[7 3.5];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print -depsc C:\D\Research\liam\SMC_EM_GLM\bernoulli\PFapprox