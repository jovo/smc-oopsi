function dgxf2(filename,fontsize,linewidth,color,fg,bg,legendfontsize,fontname,markersize,fignum)

% DGXF(FILENAME,FONTSIZE,LINEWIDTH,COLOR,FG,BG, ...
%      LEGENDFONTSIZE,FONTNAME,MARKERSIZE)
%
% Makes a nice and pretty EPS file from a matlab figure. The only
% required argument is FILENAME. If you want to set an argument
% towards the end of the list, you have to fill the arguments in
% between with empty brackets ( [] ). The arguments are as follows
% (default in parentheses):
%
% FONTSIZE (24). Size of the fonts for axis tick labels and axis
% labels.
%
% LINEWIDTH (2) Width of all of the lines in the figures.
%
% COLOR ('rgb') Color/grayscale ('rgb') or black & white figure ('bw'). (This
% essentially decides whether to use -depsc or -deps).
%
% FG ( matlab defaults ) Foreground color (default color for all
% text). Doesn't affect plot line color.
%
% BG ( matlab defaults ) Background color.
%
% LEGENDFONTSIZE (0.75*FONTSIZE) Size of the fonts for the legend.
%
% FONTNAME ('helvetica') Font name to be used. Another option is
% 'times'.
%
% MARKERSIZE (0.5*FONTSIZE) Size of plot markers.
%
% FIGNUM (gcf) Figure to operate on.
%
% Disclaimers:
% 
% 1. This script was designed to give nice EPS figures. Do
% not judge it by what it does to the figure window on the
% screen.
% 
% 2. This script will not work on figures with subplots.
%
% 3. This script is not tested on figures with titles present.
%
% By David H. Goldberg. Revised 3/11/04. 
%
% Revised by R. Jacob Vogelstein 5/8/06 to work on subplots.  Filename is
% also no longer required.  

%%% Process the input arguments

if(nargin<1)
    filename = [];
end

if(nargin<2)
  fontsize = [];
end
if(isempty(fontsize))
  fontsize=24;
end

if(nargin<3)
  linewidth = [];
end
if(isempty(linewidth))
  linewidth=2;
end

if(nargin<4)
  color = [];
end
if(isempty(color))
  color = 'rgb';
end

if(nargin<5)
  bg=[];
end

if(nargin<6)
  fg=[];
end
  
if(nargin<7)
  legendfontsize=[];
end
if(isempty(legendfontsize))
  legendfontsize=0.75*fontsize;
end

if(nargin<8)
  fontname=[];
end
if(isempty(fontname))
  fontname='helvetica';
end
 
if(nargin<9)
  markersize=[];
end
if(isempty(markersize))
  markersize=0.5*fontsize;
end

if(nargin<10)
    fignum=[];
end
if(isempty(fignum))
    fignum = gcf;
end

childaxes = get(gcf,'children');
childaxes = childaxes(find(strcmp(get(childaxes,'type'),'axes')));
for c = length(childaxes):-1:1
    
    % set the current axes to the next one in the list
    axes(childaxes(c));

    %%% Handle the axes
    axes_child_vec=get(childaxes(c),'children');

    % Adjust thickness of lines
    axes_line_idx=find(strcmp(get(axes_child_vec,'type'),'line'));
    set(axes_child_vec(axes_line_idx),'linewidth',linewidth);
    set(axes_child_vec(axes_line_idx),'markersize',markersize);

    % Adjust size of text
    axes_text_idx=find(strcmp(get(axes_child_vec,'type'),'text'));
    set(axes_child_vec(axes_text_idx),'fontsize',fontsize);
    set(axes_child_vec(axes_text_idx),'fontname',fontname);

    set(gca,'fontsize',fontsize);
    set(gca,'fontname',fontname);
    set(gca,'linewidth',linewidth);
    set(get(gca,'title'),'fontsize',fontsize);
    set(get(gca,'title'),'fontname',fontname);
    set(get(gca,'xlabel'),'fontsize',fontsize);
    set(get(gca,'xlabel'),'fontname',fontname);
    set(get(gca,'xlabel'),'verticalalignment','top');
    set(get(gca,'ylabel'),'fontsize',fontsize);
    set(get(gca,'ylabel'),'fontname',fontname);
%     set(gca,'position',[0.26 0.22 0.645 0.705]);

    %%% Do the patches

    patch_idx=find(strcmp(get(axes_child_vec,'type'),'patch'));
    set(axes_child_vec(patch_idx),'linewidth',linewidth);

    %%% If there's a legend, do that too!
    % may not work w/ multiple legends

    figure_child_vec=get(gcf,'children');
    legend_idx=find(strcmp(get(figure_child_vec,'tag'),'legend'));

    % Set the linewidth of the legend box
    set(figure_child_vec(legend_idx),'linewidth',linewidth);

    % Get the children of the legend
    legend_child_vec=get(figure_child_vec(legend_idx),'children'); 
    legend_child_line_idx=find(strcmp(get(legend_child_vec,'type'),'line'));
    set(legend_child_vec(legend_child_line_idx),'linewidth',linewidth);
    legend_child_text_idx=find(strcmp(get(legend_child_vec,'type'),'text'));
    set(legend_child_vec(legend_child_text_idx),'fontsize',legendfontsize);
    set(legend_child_vec(legend_child_text_idx),'fontname',fontname);

    %%% Set the colors
    if(not(isempty(bg)&isempty(fg)))
      set(gcf,'inverthardcopy','off');
      set(gca,'color',bg);
      set(gca,'xcolor',fg);
      set(gca,'ycolor',fg);
      set(gca,'zcolor',fg);
      set(gcf,'color',bg);
      set(figure_child_vec(legend_idx),'color',bg);
      set(figure_child_vec(legend_idx),'xcolor',fg);
      set(figure_child_vec(legend_idx),'ycolor',fg);
      set(legend_child_vec(legend_child_text_idx),'color',fg);
      set(axes_child_vec(patch_idx),'edgecolor',fg);
      set(axes_child_vec(axes_text_idx),'color',fg);
    end

    %%% Write to a postscript file
    if ~isempty(filename)
        if (strcmp(color,'bw'))
          eval(sprintf('print -deps %s',filename));
        elseif (strcmp(color,'rgb'))
          eval(sprintf('print -depsc %s',filename));
        end
    end
end