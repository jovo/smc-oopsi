function varargout = run_oopsi(F,V,P)
% this function runs our various oopsi filters, saves the results, and
% plots the inferred spike trains.  make sure that fast-oopsi and
% smc-oopsi repository are in your path if you intend to use them.  
% 
% to use the code, simply provide F, a vector of fluorescence observations,
% for each cell.  the fast-oopsi code can handle a matrix input,
% corresponding to a set of pixels, for each time bin.
%
% input
%   F: fluorescence from a single neuron
%   V: Variables necessary to define for code to function properly (optional)
%   P: Parameters of the model (optional)
%
% possible outputs
%   fast:   fast-oopsi MAP estimate of spike train, argmax_{n\geq 0} P[n|F], (fast.n) and parameter estimate (fast.P)
%   smc:    smc-oopsi estimate of {P[X_t|F]}_{t<T}, where X={n,C} or {n,C,h}, (smc.E) and parameter estimates (smc.P)

%% set code Variables

if nargin < 2,              V           = struct;   end         % create structure for algorithmic variables, if none provided
if ~isfield(V,'fast_do'),   V.fast_do  = 1;         end         % whether to use fast filter, aka, fast_oopsi
if ~isfield(V,'smc_do'),    V.smc_do   = 1;         end         % whether to use particle filter, aka, smc_oopsi
if ~isfield(V,'name'),                                          % give data a unique, time-stamped name, if there is not one specified
    cput    = clock;
    V.name  = ['oopsi_' num2str(cput(1)) '_' num2str(cput(2)) '_' num2str(cput(3)) '_' num2str(cput(4)) '_' num2str(cput(5))];
end         
V.name_dat = ['~/Research/oopsi/smc-oopsi/data/jovo/' V.name];  % filename for data
V.name_fig = ['~/Research/oopsi/smc-oopsi/figs/jovo/' V.name];  % filename for figure
if nargin < 3,              P           = struct;   end         % create structure for parameters, if none provided
save(V.name_dat,'V')

%% infer spikes and estimate parameters

if V.fast_do == 1                                          % infer spikes using fast-oopsi
    [fast.n fast.P fast.V]= fast_oopsi(F,V,P);                                
    save(V.name_dat,'fast','-append')
end

if V.smc_do == 1                                           % infer spikes using smc-oopsi
    siz=size(F); if siz(1)>1, F=F'; end
    if V.fast_do == 1; 
        V.fast_n= fast.n;
        P.alpha = fast.P.a;
        P.beta  = fast.P.b;
%         P.k     = f(P.lam);
        P.gamma = 0;
        P.zeta  = fast.P.sig;
        P.tau_c = V.dt/(1-fast.P.gam);
    end
    [smc.E smc.P smc.V] = smc_oopsi(F,V,P);
    save(V.name_dat,'smc','-append');
end

%% provide outputs for later analysis

if nargout == 1
    if V.fast_do == 1; 
        varargout{1} = fast;
    else
        varargout{1} = smc;
    end
else
    varargout{1} = fast;
    varargout{2} = smc;
end

%% plot results

fig=figure(3); clf,
V.T     = length(F);
if V.fast_do==1, V.dt=fast.V.dt; else  V.dt=smc.V.dt; end
if (V.fast_do==1 && V.smc_do==1), nrows=3; else nrows=2; end
gray    = [.75 .75 .75];            % define gray color
col   = [1 0 0; 0.2 0.2 1];         % define colors for mean
ccol  = col+.4; ccol(ccol>1)=1;     % define colors for std
inter   = 'tex';                    % interpreter for axis labels
xlims   = [45 V.T-2];               % xmin and xmax for current plot
fs      = 12;                       % font size
ms      = 5;                        % marker size for real spike
sw      = 2;                        % spike width
lw      = 2;                        % line width
xticks  = 0:1/V.dt:V.T;             % XTick positions
skip    = round(length(xticks)/5);
xticks  = xticks(1:skip:end);
tvec_o=xlims(1):xlims(2);
i=0;
if isfield(V,'true_n'), spt = find(V.true_n==1); end               % find spike times

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
if V.fast_do==1
    i=i+1; h(i)=subplot(nrows,1,i); hold on,
    n_fast=fast.n/max(fast.n);
    spts=find(n_fast>1e-3);
    stem(spts,n_fast(spts),'Marker','none','LineWidth',sw,'Color','k')
    if isfield(V,'true_n'),
        stem(spt,V.true_n(spt)+0.1,'Marker','v','MarkerSize',ms,'LineStyle','none','Color',gray,'MarkerFaceColor',gray,'MarkerEdgeColor',gray)
    end
    axis([xlims 0 1.1])
    hold off,
    ylab=ylabel([{'Fast'}; {'Filter'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    box off
end

% plot smc-oopsi output
if V.smc_do == 1
    i=i+1; h(i)=subplot(nrows,1,i); hold on,
    spts=find(smc.E.nbar>1e-3);
    stem(spts,smc.E.nbar(spts),'Marker','none','LineWidth',sw,'Color','k')
    if isfield(V,'true_n'),
        stem(spt,V.true_n(spt)+0.1,'Marker','v','MarkerSize',ms,'LineStyle','none','Color',gray,'MarkerFaceColor',gray,'MarkerEdgeColor',gray)
    end
    axis([xlims 0 1.1])
    hold off,
    ylab=ylabel([{'Particle'}; {'Filter'}],'Interpreter',inter,'FontSize',fs);
    set(ylab,'Rotation',0,'HorizontalAlignment','right','verticalalignment','middle')
    set(gca,'YTick',0:2,'YTickLabel',[])
    set(gca,'XTick',xticks,'XTickLabel',[])
    box off
end

% label last subplot
set(gca,'XTick',xticks,'XTickLabel',round(xticks*V.dt*100)/100)
xlabel('Time (sec)','FontSize',fs)
linkaxes(h,'x')

% print fig
wh=[7 3];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
print('-depsc',V.name_fig)
print('-dpdf',V.name_fig)
saveas(fig,V.name_fig)
