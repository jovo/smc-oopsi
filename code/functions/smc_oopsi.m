function [S P] = smc_oopsi(F,V,P0)
% this function prepares things for running our smc-oopsi code. simply
% provide F, a vector of fluorescence observations (unfiltered), and the
% code will output a best guess of the spike train.  parameters are
% initialized (by default) using the fast_oopsi code, so that function must
% be in the path as well.  
% 
% input
%   F: 1D fluorescence trace from a single neuron
%   V: structure of variables necessary to define for code to function properly (optional)
%   P0: model parameters (optional)
% 
% output
%   S: smc particle weights and positions
%   P: parameter estimates

%% set code Variables

if nargin < 2, V = struct; end
if ~isfield(V,'T'),      V.T = length(F); end           % # of time steps
if ~isfield(V,'freq'),   V.freq = 1; end                % # time steps between observations
if ~isfield(V,'T_o'),    V.T_o = V.T; end               % # of observations
if ~isfield(V,'N'),      V.N = 99; end                  % # particles
if ~isfield(V,'M'),      V.M = 0; end                   % # of spike history terms
if ~isfield(V,'pf'),     V.pf = 1; end                  % whether to use conditional sampler
if ~isfield(V,'x'),      V.x = ones(1,V.T); end         % stimulus
if ~isfield(V,'name'),   V.name='oopsi_output'; end     % name for output and figure
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
    if ~isfield(V,'C_params'), V.C_params   = 1; end    % calcium params
    if ~isfield(V,'n_params'), V.n_params   = 1; end    % b,k
    if ~isfield(V,'h_params'), V.h_params   = 1; end    % w
    if ~isfield(V,'F_params'), V.F_params   = 1; end    % alpha, beta
    if ~isfield(V,'G_params'), V.G_params   = 1; end    % gamma, zeta
    if ~isfield(V,'Scan'),     V.Scan       = 0; end    % epi or scan
    if ~isfield(V,'FastInit'), V.FastInit   = 1; end    % epi or scan
    if V.FastInit == 1; if ~isfield(V,'ptile'), V.ptile = 0.95; end; end
end

% % V.TrueSpk = double(n);
% V.holdTau=1;

%% initialize model Parameters

if nargin < 3, P0 = struct; end
if ~isfield(P0,'k'),     P0.k     = log(-log(1-100/V.T)/V.dt); end  % linear filter
if ~isfield(P0,'tau_c'), P0.tau_c = 1; end                          % calcium decay time constant (sec)
if ~isfield(P0,'A'),     P0.A     = 50; end                         % change ins [Ca++] after a spike (\mu M)
if ~isfield(P0,'C_0'),   P0.C_0   = 0; end                          % baseline [Ca++] (\mu M)
if ~isfield(P0,'C_init'),P0.C_init= 0; end                          % initial [Ca++] (\mu M)
if ~isfield(P0,'sigma_c'),P0.sigma_c= 0.1; end                        % standard deviation of noise (\mu M)
if ~isfield(P0,'n'),     P0.n     = 1; end                          % hill equation exponent
if ~isfield(P0,'k_d'),   P0.k_d   = 200; end                       % hill coefficient
if ~isfield(P0,'alpha'), P0.alpha = mean(F); end                    % scale of F
if ~isfield(P0,'beta'),  P0.beta  = min(F); end                     % offset of F
if ~isfield(P0,'gamma'), P0.gamma = 0; end                          % scaled variance
if ~isfield(P0,'zeta'),  P0.zeta  = P0.alpha/5; end                 % constant variance
if ~isfield(P0,'a'),     P0.a     = V.dt/P0.tau_c; end          

if V.M==1                                                           % if there are spike history terms
    if ~isfield(P0,'omega'),   P0.omega   = -1; end                 % weight
    if ~isfield(P0,'tau_h'),   P0.tau_h   = 0.02; end               % time constant
    if ~isfield(P0,'sigma_h'), P0.sigma_h = 0; end                  % stan dev of noise
end

%% infer spikes and estimate parameters

siz=size(F); if siz(1)>1, F=F'; end
[S P] = GOOPSI_main_v1_0(F,P0,V);
save(V.name,'S','P','F','V');

%% plot results

fig     = figure(3); clf,
if V.FastInit==0; nrows=2; else nrows=3; end
gray    = [.75 .75 .75];            % define gray color
col   = [1 0 0; 0.2 0.2 1];         % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter   = 'tex';                    % interpreter for axis labels
xlims   = [45 V.T-2*V.freq];    % xmin and xmax for current plot
fs      = 12;                       % font size
ms      = 5;                       % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
xticks  = 0:1/V.dt:V.T;               % XTick positions
skip    = round(length(xticks)/5);
xticks  = xticks(1:skip:end);
tvec_o=xlims(1):V.freq:xlims(2);
i=0;
if isfield(V,'n'), spt = find(V.n==1); end               % find spike times

% plot fluorescence data
i=i+1; h(i)=subplot(nrows,1,i); hold on
plot(tvec_o,z1(F(tvec_o)),'.-k','LineWidth',2,'MarkerSize',ms*.75);
% stem(V.n,'Marker','none','LineWidth',sw,'Color','k')
ylab=ylabel([{'Fluorescence'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',[],'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
axis([xlims 0 1.1])


% plot fast-oopsi output 
if V.FastInit==1
    i=i+1; h(i)=subplot(nrows,1,i); hold on,
    n_fast=S.n_fast/max(S.n_fast);
    spts=find(n_fast>1e-3);
    stem(spts,n_fast(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
    if isfield(V,'n'),
        stem(spt,V.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
    end
    axis([xlims 0 1])
    hold off,
    ylab=ylabel([{'Fast'}; {'Filter'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    box off
end

% plot smc-oopsi output 
i=i+1; h(i)=subplot(nrows,1,i); hold on,
spts=find(S.nbar>1e-3);
stem(spts,S.nbar(spts),'Marker','none','LineWidth',sw,'Color',col(2,:))
if isfield(V,'n'), 
    stem(spt,V.n(spt),'Marker','v','MarkerSize',ms,'LineStyle','none','Color','k','MarkerFaceColor','k','MarkerEdgeColor','k')
end
axis([xlims 0 1])
hold off,
ylab=ylabel([{'Particle'}; {'Filter'}],'Interpreter',inter,'FontSize',fs);
set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
set(gca,'YTick',0:2,'YTickLabel',[])
set(gca,'XTick',xticks,'XTickLabel',[])
box off

% label last subplot
set(gca,'XTick',xticks,'XTickLabel',round(xticks*V.dt*100)/100)
xlabel('Time (sec)','FontSize',fs)
linkaxes(h,'x')

% print fig
wh=[7 3];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-depsc',V.name)
print('-dpdf',V.name)
