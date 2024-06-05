function plotGivenFieldVariables(folderName, scaleFactor, varName, step, gridSize)
a=load([folderName,'\pd_model.mat']);
coor = a.pd_model.coor;

resFile = [folderName,'\result_step_',num2str(step),'.mat'];
if(exist(resFile,'file'))
    b = load(resFile);
end
% u1 = b.result.u1;
% u2 = b.result.u2;
u1 = 0;
u2 = 0;
eval(['val=b.result.',varName,';']);
clf
scatter(coor(:,1)+scaleFactor*u1, coor(:,2)+scaleFactor*u2,gridSize,val,'filled')
colormap(jet)
axis equal
axis off
colorbar
if(strcmp(varName,'dmg'))
    caxis([0,1])
end
title(varName)
set(gca,'fontsize',16,'fontname','times new roman')
end