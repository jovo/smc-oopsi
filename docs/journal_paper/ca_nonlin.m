function ca_nonlin3(o)
% F ~ N[f(C),g(C)]

[n k_d alpha beta gamma zeta] = GetParams;

if o > alpha+beta || o < beta
    error('that''s not possible')
end

% compute mean
finv    = ((k_d*(o-beta))./(alpha-o+beta)).^(1/n); %initialize search with f^{-1}(o)
options = optimset('GradObj','off','Hessian','off');
ffit    = fminunc(@fnlogL,finv,options);     %max P(O|H)

% compute variance
syms Cest
logL=-fnlogL(Cest);
dlogL=diff(logL,'Cest');            %dlog L / dC
ddlogL=diff(dlogL,'Cest');          %ddlog L / dCC
VC = -1/ddlogL;                     %neg inverse
Cest = ffit;                        %eval at max P(O|H)
varF=eval(VC);                      %variance approximation

    function mu_F = fmu_F(C)        %this function compute E[F]=f(C)
        mu_F    = alpha*C.^n./(C.^n+k_d)+beta; 
    end

    function var_F = fvar_F(C)      %this function compute V[F]=f(C)
        var_F   = gamma*C.^n./(C.^n+k_d)+zeta;
    end

    function [logL dlogL ddlogL] = fnlogL(C)        %this function compute log L = log P(O|H)
        logL = (((o-fmu_F(C)).^2)./(2*fvar_F(C))+log(fvar_F(C)))/2;
        if nargout > 1
                dlogL=-(o-alpha*C^n/(C^n+k_d)-beta)/(2*gamma*C^n/(C^n+k_d)+2*zeta)*(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)+1/2*(o-alpha*C^n/(C^n+k_d)-beta)^2/(2*gamma*C^n/(C^n+k_d)+2*zeta)^2*(2*gamma*C^n*n/C/(C^n+k_d)-2*gamma*(C^n)^2/(C^n+k_d)^2*n/C)-1/2*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)/(gamma*C^n/(C^n+k_d)+zeta);
            if nargout > 2
                ddlogL=-(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)^2/(2*gamma*C^n/(C^n+k_d)+2*zeta)+2*(o-alpha*C^n/(C^n+k_d)-beta)/(2*gamma*C^n/(C^n+k_d)+2*zeta)^2*(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)*(2*gamma*C^n*n/C/(C^n+k_d)-2*gamma*(C^n)^2/(C^n+k_d)^2*n/C)-(o-alpha*C^n/(C^n+k_d)-beta)/(2*gamma*C^n/(C^n+k_d)+2*zeta)*(-alpha*C^n*n^2/C^2/(C^n+k_d)+alpha*C^n*n/C^2/(C^n+k_d)+3*alpha*(C^n)^2*n^2/C^2/(C^n+k_d)^2-2*alpha*(C^n)^3/(C^n+k_d)^3*n^2/C^2-alpha*(C^n)^2/(C^n+k_d)^2*n/C^2)-(o-alpha*C^n/(C^n+k_d)-beta)^2/(2*gamma*C^n/(C^n+k_d)+2*zeta)^3*(2*gamma*C^n*n/C/(C^n+k_d)-2*gamma*(C^n)^2/(C^n+k_d)^2*n/C)^2+1/2*(o-alpha*C^n/(C^n+k_d)-beta)^2/(2*gamma*C^n/(C^n+k_d)+2*zeta)^2*(2*gamma*C^n*n^2/C^2/(C^n+k_d)-2*gamma*C^n*n/C^2/(C^n+k_d)-6*gamma*(C^n)^2*n^2/C^2/(C^n+k_d)^2+4*gamma*(C^n)^3/(C^n+k_d)^3*n^2/C^2+2*gamma*(C^n)^2/(C^n+k_d)^2*n/C^2)-1/2*(gamma*C^n*n^2/C^2/(C^n+k_d)-gamma*C^n*n/C^2/(C^n+k_d)-3*gamma*(C^n)^2*n^2/C^2/(C^n+k_d)^2+2*gamma*(C^n)^3/(C^n+k_d)^3*n^2/C^2+gamma*(C^n)^2/(C^n+k_d)^2*n/C^2)/(gamma*C^n/(C^n+k_d)+zeta)+1/2*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)^2/(gamma*C^n/(C^n+k_d)+zeta)^2;            end
        end
    end

tempo = o;
syms Cest alpha beta n k_d gamma zeta o
logL=-fnlogL(Cest);
dlogL=diff(logL,'Cest');            %dlog L / dC
ddlogL=diff(dlogL,'Cest');          %ddlog L / dCC

[n k_d alpha beta gamma zeta] = GetParams;
C       = ffit;                     %eval at max P(O|H)
o       = tempo;
varF2=-1/(-(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)^2/(2*gamma*C^n/(C^n+k_d)+2*zeta)+2*(o-alpha*C^n/(C^n+k_d)-beta)/(2*gamma*C^n/(C^n+k_d)+2*zeta)^2*(-alpha*C^n*n/C/(C^n+k_d)+alpha*(C^n)^2/(C^n+k_d)^2*n/C)*(2*gamma*C^n*n/C/(C^n+k_d)-2*gamma*(C^n)^2/(C^n+k_d)^2*n/C)-(o-alpha*C^n/(C^n+k_d)-beta)/(2*gamma*C^n/(C^n+k_d)+2*zeta)*(-alpha*C^n*n^2/C^2/(C^n+k_d)+alpha*C^n*n/C^2/(C^n+k_d)+3*alpha*(C^n)^2*n^2/C^2/(C^n+k_d)^2-2*alpha*(C^n)^3/(C^n+k_d)^3*n^2/C^2-alpha*(C^n)^2/(C^n+k_d)^2*n/C^2)-(o-alpha*C^n/(C^n+k_d)-beta)^2/(2*gamma*C^n/(C^n+k_d)+2*zeta)^3*(2*gamma*C^n*n/C/(C^n+k_d)-2*gamma*(C^n)^2/(C^n+k_d)^2*n/C)^2+1/2*(o-alpha*C^n/(C^n+k_d)-beta)^2/(2*gamma*C^n/(C^n+k_d)+2*zeta)^2*(2*gamma*C^n*n^2/C^2/(C^n+k_d)-2*gamma*C^n*n/C^2/(C^n+k_d)-6*gamma*(C^n)^2*n^2/C^2/(C^n+k_d)^2+4*gamma*(C^n)^3/(C^n+k_d)^3*n^2/C^2+2*gamma*(C^n)^2/(C^n+k_d)^2*n/C^2)-1/2*(gamma*C^n*n^2/C^2/(C^n+k_d)-gamma*C^n*n/C^2/(C^n+k_d)-3*gamma*(C^n)^2*n^2/C^2/(C^n+k_d)^2+2*gamma*(C^n)^3/(C^n+k_d)^3*n^2/C^2+gamma*(C^n)^2/(C^n+k_d)^2*n/C^2)/(gamma*C^n/(C^n+k_d)+zeta)+1/2*(gamma*C^n*n/C/(C^n+k_d)-gamma*(C^n)^2/(C^n+k_d)^2*n/C)^2/(gamma*C^n/(C^n+k_d)+zeta)^2);

% just find max to check numerically
C       = 0.01:.001:50;             %possible values for C
mu_F    = fmu_F(C);                 %E[F] forall C
var_F   = fvar_F(C);                %V[F] forall C
L       = -fnlogL(C);                %L eval at o forall C
[foo ind]=max(L);                   %get max numerically to check
fnum    = C(ind);                   

% check approx quality
Lapprox = 1/sqrt(2*pi*varF)*exp(-.5*(ffit-C).^2/varF);
Lapprox2 = 1/sqrt(2*pi*varF2)*exp(-.5*(ffit-C).^2/varF2);

%% make plots
yfs = 14;
xfs = 14;
ms  = 25;

% plot E[F] vs. log C
figure(11); clf, 
subplot(311)
semilogx(C,mu_F,'k','linewidth',2);  grid on; hold on
finv=((k_d*(o-beta))./(alpha-o+beta)).^(1/n);
semilogx(finv,o,'.k','linewidth',2,'MarkerSize',ms)

axis([min(C) max(C) 0 1]);  
set(gca,'fontsize',xfs)
set(gca,'XTickLabel',[])
ylab=['mu_v'];
set(get(gca,'YLabel'),'String',texlabel(ylab),'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs)
% get(gca,'YTick')
% set(gca,'YTick',[0:5:1])
% xlab=xlabel('[Ca^{2+}]');
% set(xlab,'fontsize',yfs)

% plot var[F] vs. log C
subplot(312);
semilogx(C,(var_F),'k','linewidth',2); grid on; hold on
varf=gamma*finv.^n./(finv.^n+k_d)+zeta;
semilogx(finv,varf,'.k','linewidth',2,'MarkerSize',ms)

axis([min(C) max(C) 0.00100 0.00110])
set(gca,'fontsize',xfs)
ylab=['sigma_v^2'];
set(get(gca,'YLabel'),'String',texlabel(ylab),'fontsize',yfs,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')

% axis('tight')
% YTick=get(gca,'YTick');
% set(gca,'YTick',[0.00102 0.00108])
% xlab=xlabel('[Ca^{2+}]');
% set(xlab,'fontsize',yfs)
% ylab=ylabel('F');
% xlab=['sigma_{Ca^{2+}}^2'];
% set(get(gca,'XLabel'),'String',texlabel(xlab),'fontsize',yfs)
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle','fontsize',yfs)

% plot L vs. C
subplot(313); hold on; %title(['o=' num2str(o) ', mu=' num2str(ffit) ', var=' num2str(varF2)])
plot(C,(exp(L)+o)/max(exp(L)+o),'k','linewidth',2); 
ylab=ylabel({'Normalized';'Likelihood'});
% set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(ylab,'fontsize',yfs)
xlab=xlabel('[Ca^2^+]');
set(xlab,'fontsize',yfs)
plot(C, Lapprox/max(Lapprox),'--k','linewidth',2); 
set(gca,'fontsize',xfs)
% plot(C, Lapprox2/max(Lapprox2),'-.r'); 
axis([ffit-3*sqrt(varF2) ffit+3*sqrt(varF2) 0 1.02])

% make eps
fig=figure(11);
wh=2*[3 3];
set(fig,'PaperPosition',[0 11-wh(2) wh]);
print('-deps','ca_nonlin')

end

function [n k_d alpha beta gamma zeta] = GetParams % set parameters

n       = 1.2;                      %hill coefficient
k_d     = 1.3;                      %dissociation constant
alpha   = 1;                        %mean gain
beta    = 0;                       %mean offset
gamma   = 1e-4;             %var gainh
zeta    = 1e-3;            %var offset

end