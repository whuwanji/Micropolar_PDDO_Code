function makeGifU(folderName,fileNameOut)
a=load([folderName,'\pd_model.mat']);
coor = a.pd_model.coor;
nt = a.pd_model.nt;
filename=[fileNameOut,'.gif'];
for i = 100:10:3000
    resFile = [folderName,'\result_step_',num2str(i),'.mat'];
    if(exist(resFile,'file'))
        b = load(resFile);
    else
        break;
    end
   
    u1 = b.result.u1;
    u2 = b.result.u2;
    Fig = figure(1);
    scatter(coor(:,1)+5*u1,coor(:,2)+5*u2,10,u2,'filled')
    colormap(jet)
    axis equal
    axis off
    colorbar
    %caxis([0,1])
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
    title('u_2')
    set(gca,'fontsize',16,'fontname','times new roman')
    pause(0.00001)
    frame = getframe(Fig); 

    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256);
    if i == 100
         imwrite(imind,cm,filename,'gif','WriteMode','overwrite', 'Loopcount',inf);
    else
         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.02);
   end
end
