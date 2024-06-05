function makeGif(folderName,fileNameOut, scaleFactor, varName, startStep, step, endStep, gridSize)
a=load([folderName,'\pd_model.mat']);
coor = a.pd_model.coor;
nt = a.pd_model.nt;
filename=[fileNameOut,'.gif'];
for i = startStep:step:endStep
    resFile = [folderName,'\result_step_',num2str(i),'.mat'];
    if(exist(resFile,'file'))
        b = load(resFile);
    else
        break;
    end
    u1 = b.result.u1;
    u2 = b.result.u2;
    val = b.result.dmg;
    eval(['val=b.result.',varName,';']);
    Fig = figure(1);
    clf
    scatter(coor(:,1)+scaleFactor*u1, coor(:,2)+scaleFactor*u2,gridSize,val,'filled')
    colormap(jet)
    axis equal
    axis off
    colorbar
    if(strcmp(varName,'dmg'))
        caxis([0,1])
    end
    %     xlabel('x')
    %     ylabel('y')
    %     zlabel('z')
    title(varName)
    set(gca,'fontsize',16,'fontname','times new roman')
    pause(0.00001)
    frame = getframe(Fig);
    
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == startStep
        imwrite(imind,cm,filename,'gif','WriteMode','overwrite', 'Loopcount',inf);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.02);
    end
end
