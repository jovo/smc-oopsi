function [name, file]=myfig(gcf,name,file);
% [name, file]=myfig(gcf,name,file);
% take input gcf and set paperpositionmode to auto and all tickmodes to manual
% name should include directory 
% print an eps file and a fig file
% additionally copy the m.file used to make the figure there

in=input(['do you really want to save the figure as '...
	,name,'.eps? [y for yes, else no]'],'s');

if strcmp(in,'y')

	set(gcf,'paperpositionmode','auto')
	t=get(gcf,'children');
	for k=1:length(t);
		set(t(k),'xtickmode','manual','ytickmode','manual','ztickmode','manual');
	end

	print(gcf,'-depsc',strcat(name,'.eps'))
	saveas(gcf,strcat(name,'.fig'));
	%save(strcat(name,'.mat'))

	eval(['!gv ',name,'.eps']);
	if exist('file'); eval(['!cp ',file,' ',name,'.m']);end

end