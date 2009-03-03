% this script does makes a figure showing how the smc-em approach
% can do superresolution in the temporal domain
% 
% 1) set simulation metadata (eg, dt, T, # particles, etc.)
% 2) initialize parameters
% 3) generate fake data for an array of subsampling and noisiness
% 4) infers spikes using a variety of approaches
% 5) plots results

clear, clc, close all, fprintf('\nArray Simulation Fig\n')

%% 1) set simulation metadata

Sim.T       = 50;                                   % # of time steps
Sim.dt      = 1/50;                                 % time step size
Sim.freq    = 5;                                    % # of time steps between observations
Sim.Nsec    = Sim.T*Sim.dt;                         % # of actual seconds
Sim.T_o     = Sim.T;                                % # of observations
Sim.tvec    = Sim.dt:Sim.dt:Sim.Nsec;               % vector of times
Sim.N       = 200;                                  % # of particles
Sim.M       = 0;                                    % number of spike history dimensions
Sim.pf      = 1;                                    % use conditional sampler (not prior) when possible
Sim.StimDim = 1;                                    % # of stimulus dimensions
Sim.x       = ones(1,Sim.T);                        % stimulus

Sim.Mstep   = false;                                % do M-step
Sim.C_params = true;                                 % whether to estimate calcium parameters {tau,A,C_0,sig}
Sim.n_params = true;                                 % whether to estimate rate governing parameters {b,k}
Sim.h_params = false;                                % whether to estimate spike history parameters {h}
Sim.F_params = false;                                % whether to estimate observation parameters {alpha,beta,gamma,zeta}
Sim.MaxIter = 0;                                    % max # of EM iterartions

%% 2) initialize parameters

% initialize barrier and wiener filter parameters
spt     = [12 27 38 43];                            % spike times
P.rate  = numel(spt)/(Sim.T*Sim.dt);                % expected spike rate
P.A     = 1;                                        % jump size ($\mu$M)
P.tau   = 2;                                        % calcium decay time constant (sec)
P.lam   = Sim.T/(P.rate*P.A)*Sim.dt;                % expected jump size ber time bin
P.sig   = 1;                                        % standard deviation of noise (\mu M)

% initialize particle filter parameters
P.k         = log(-log(1-P.rate*Sim.dt)/Sim.dt);    % linear filter
P.tau_c     = P.tau;
P.A         = 10;
P.C_0       = 0.1;                                   % baseline [Ca++]
P.C_init    = P.C_0;                                % initial [Ca++]
P.sigma_c   = P.sig;
P.n         = 1.0;                                  % hill equation exponent
P.k_d       = 200;                                  % hill coefficient
P.alpha     = 1;                                    % F_max
P.beta      = 0;                                    % F_min
P.a         = Sim.dt/P.tau_c;
gamma       = 1e-5;                                 % scaled variance
zeta        = 4*gamma;                              % constant variance

%% 3) simulate data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% n=rand(Sim.T,1)<1-exp(-exp(P.k)*Sim.dt);
n=zeros(Sim.T,1); n(spt)=1;
Sim.n = double(n); Sim.n(Sim.n==0)=NaN;             % for plotting purposes in ParticleFiltD
C=zeros(Sim.T,1);
epsilon_c = P.sigma_c*sqrt(Sim.dt)*randn(1,Sim.T);      %generate noise on calcium
for t=2:Sim.T                                           %recursively update calcium
    C(t)  = (1-Sim.dt/P.tau_c)*C(t-1) + P.a*P.C_0+ P.A*n(t) + epsilon_c(t);   
end
S=Hill_v1(P,C);
vec=[1 2 4 8];
Algs=7; Sim.Alg = 7;

for sub=1:length(vec)
    Sim.freq    = vec(sub);
    Sim.T_o     = Sim.T/Sim.freq;                       %fix # of observations
    obs         = Sim.freq:Sim.freq:Sim.T;
    for nois=1:length(vec)
        noi     = vec(nois);
        P.gamma = noi*gamma;
        P.zeta  = noi*zeta;
        F_temp  = P.alpha*S+P.beta+sqrt(P.gamma*S+P.zeta).*randn(Sim.T,1);
        F       = zeros(size(S));
        F(obs)  = F_temp(obs);
        F(F<=0) = eps;
        figure, plot(F+1,'.-'); hold on, stem(n)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 4) infer spikes and estimate parameters
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fprintf('\nfreq=%g, noise=%g\n',Sim.freq,noi)
        I{sub,nois} = DataComp13(F,P,Sim);
        I{sub,nois}.F = F;
    end
end
save ArraySim

%% 5) plot results

load('ArraySim.mat')

fig=figure(100); clf,
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
        stem(spts,BarVar(spts),'Marker','none','LineWidth',sw,'Color',ccol(2,:));
        spts=find(I{sub,nois}.M.nbar>1e-3);
        stem(spts,I{sub,nois}.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
        stem(spt,Sim.n(spt)-.05,'Marker','v','MarkerSize',ms2,'LineStyle','none','MarkerFaceColor','k','MarkerEdgeColor','k')
        axis([xlims 0 2+yshift])
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
print('-depsc','ArraySim')

%% 6) plot b&w results

load('ArraySim.mat')

fig=figure(100); clf,
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
yshift = .15;

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