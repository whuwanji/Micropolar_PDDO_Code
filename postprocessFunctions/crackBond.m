function crackLen = crackBond(folderName, steps)
a=load([folderName,'\pd_model.mat']);
cLen = a.pd_model.crackLen;
coor = a.pd_model.coor;
nt = numel(steps);
crackLen = NaN(nt,2);
for i = 1:1:nt
    resFile = [folderName,'\result_step_',num2str(steps(i)),'.mat'];
    if(exist(resFile,'file'))
        b = load(resFile);
    else
        break;
    end
    
    dmg = b.result.dmg;
%         figure(1)
%         clf
%         scatter(coor(:,1),coor(:,2),5,dmg,'filled')
%         colormap(jet)
%         axis equal
%         axis off
%         pause(0.001)
    crackLen(i,1) =  b.result.apply_disp_u2;
    if(i==1)
        crackLen0 = max(coor(dmg>0.1,1));
    end  
    crackLen(i,2) = max(coor(dmg>0.1,1)) - crackLen0;
end