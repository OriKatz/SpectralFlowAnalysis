function [] = MySaveFig(f,folder,filename)


% set(gcf, 'Position', get(0, 'Screensize'));
set(gcf, 'Position',[1,1,1920,1080]);
% set(gcf, 'Position',[1,1,1700,500]);

savefig(folder+filename+".fig");
saveas(f,folder+filename+".eps",'epsc');
saveas(f,folder+filename+".jpg");

end

