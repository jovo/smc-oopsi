function [S P] = smc_oopsi(F,V,P0)
% this function prepares things for running our smc-oopsi code
% 
% input
%   F: 1D fluorescence trace from a single neuron
%   V: structure of variables necessary to define for code to function properly
%   P0: model parameters
% 
% output
%   S: smc particle weights and positions
%   P: parameter estimates

%% set code Variables

if nargin < 2, V = struct; end
if ~isfield(V,'T'),      V.T = length(F); end           % # of time steps
if ~isfield(V,'freq'),   V.freq = 1; end                % # time steps between observations
if ~isfield(V,'T_o'),    V.T_o = V.T; end               % # of observations
if ~isfield(V,'N'),      V.N = 20; end                  % # particles
if ~isfield(V,'M'),      V.M = 0; end                   % # of spike history terms
if ~isfield(V,'pf'),     V.pf = 1; end                  % whether to use conditional sampler
if ~isfield(V,'x'),      V.x = ones(1,V.T); end         % stimulus
if ~isfield(V,'dt'),                                    % frame duration
    fr = input('what was the frame rate for this movie (in Hz)? ');
    V.dt = 1/fr;
end
if ~isfield(V,'MaxIter'),                               % max # of iterations to estimate params
    reply = input('do you want to estimate parameters? y/n [y] (case sensitive): ', 's');
    if reply == 'y'; V.MaxIter = 10;
    else V.MaxIter = 0; end
end

% set which parameters to estimate
if V.MaxIter>1;
    if ~isfield(V,'C_params'), V.C_params   = 0; end    % calcium params
    if ~isfield(V,'n_params'), V.n_params   = 1; end    % b,k
    if ~isfield(V,'h_params'), V.h_params   = 1; end    % w
    if ~isfield(V,'F_params'), V.F_params   = 1; end    % alpha, beta
    if ~isfield(V,'G_params'), V.G_params   = 1; end    % gamma, zeta
    if ~isfield(V,'Scan'),     V.Scan       = 0; end    % epi or scan
    if ~isfield(V,'FastInit'), V.FastInit   = 1; end    % epi or scan
end

% % V.TrueSpk = double(n);
% V.holdTau=1;

%% initialize model Parameters

if nargin < 3, P0 = struct; end
if ~isfield(P0,'k'),     P0.k     = log(-log(1-V.dt)/V.dt); end   % linear filter
if ~isfield(P0,'tau_c'), P0.tau_c = 1; end                        % calcium decay time constant (sec)
if ~isfield(P0,'A'),     P0.A     = 50; end                       % change ins [Ca++] after a spike (\mu M)
if ~isfield(P0,'C_0'),   P0.C_0   = 0; end                        % baseline [Ca++] (\mu M)
if ~isfield(P0,'C_init'),P0.C_init= 0; end                        % initial [Ca++] (\mu M)
if ~isfield(P0,'sigma_c'),P0.sigma_c= 0; end                      % standard deviation of noise (\mu M)
if ~isfield(P0,'n'),     P0.n     = 1; end                        % hill equation exponent
if ~isfield(P0,'k_d'),   P0.k_d   = 1000; end                     % hill coefficient
if ~isfield(P0,'alpha'), P0.alpha = mean(F); end                        % scale of F
if ~isfield(P0,'beta'),  P0.beta  = min(F); end                        % offset of F
if ~isfield(P0,'gamma'), P0.gamma = 0; end                        % scaled variance
if ~isfield(P0,'zeta'),  P0.zeta  = P0.alpha/5; end                % constant variance
if ~isfield(P0,'a'),     P0.a     = V.dt/P0.tau_c; end          

if V.M==1                                                     % if there are spike history terms
    if ~isfield(P0,'omega'),   P0.omega   = -1; end               % weight
    if ~isfield(P0,'tau_h'),   P0.tau_h   = 0.02; end             % time constant
    if ~isfield(P0,'sigma_h'), P0.sigma_h = 0; end                % stan dev of noise
end



%% 4) infer spikes and estimate parameters

[S P] = GOOPSI_main_v1_0(F,P0,V);

%% 5) plot results

fig     = figure(2); clf,
nrows=4;
if V.M>0, nrows=nrows+1; end
if var(V.x)~=0, nrows=nrows+2; end
gray    = [.75 .75 .75];            % define gray color
col   = [1 0 0; 0.2 0.2 1];         % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter   = 'tex';                    % interpreter for axis labels
xlims   = [45 V.T-2*V.freq];    % xmin and xmax for current plot
fs      = 12;                       % font size
ms      = 20;                       % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
xticks  = 0:round(V.T/5):V.T;               % XTick positions
spt = find(V.n==1);               % find spike times
tvec_o=xlims(1):V.freq:xlims(2);
i=0;

% plot stimulus
if var(V.x)~=0
    i=i+1; subplot(nrows,1,i), hold on
    plot(V.x,'k','LineWidth',2);
    ylab=ylabel([{'External'}; {'Stimulus'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[],'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    axis([xlims min(V.x*.9) max(V.x*1.1)])
end

% plot spike history
if V.M>0
    i=i+1; subplot(nrows,1,i), hold on
    plot(P0.omega*h,'k','LineWidth',2);
    ylab=ylabel([{'Simulated'}; {'Spike History'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[],'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    axis([xlims min(P0.omega*h) 0.1])
end

% plot prob spiking
if var(V.x)~=0
    i=i+1; subplot(nrows,1,i), hold on
    plot(p,'k','LineWidth',2);
    ylab=ylabel([{'Probability'}; {'of Spiking'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',[],'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    axis([xlims 0 1])
end

% plot spike train
i=i+1; subplot(nrows,1,i), hold on
stem(spt,V.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
ylab=ylabel([{'Simulated'}; {'Spike Train'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1])

% plot calcium
i=i+1; subplot(nrows,1,i), hold on
plot(C/P0.k_d,'k','LineWidth',2);
ylab=ylabel([{'Simulated'}; {'Calcium'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1*max(C/P0.k_d)])

% plot fluorescence data
i=i+1; subplot(nrows,1,i), hold on
plot(tvec_o,z1(F(tvec_o)),'.-k','LineWidth',2,'MarkerSize',ms*.75);
% stem(V.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Simulated'}; {'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])

% plot inferred spike trains
i=i+1; subplot(nrows,1,i), hold on,
stem(spt,V.n(spt),'Marker','.','MarkerSize',ms,'LineWidth',sw,'Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
spts=find(I.M.nbar>1e-3);
stem(spts,I.M.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
axis([xlims 0 1])
hold off,
ylab=ylabel([{'SMC'}; {'OOPSI'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
box off

% label last subplot
set(gca,'XTick',xticks,'XTickLabel',(xticks-xlims(1))*V.dt)
xlabel('Time (sec)','FontSize',fs)