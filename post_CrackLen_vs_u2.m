clc; clear
addpath('postprocessFunctions')
figure(1); clf
folderName = {'resFolder'};
lcArr = 0.5;
color = 'rgb';
shape = {'-','--','-.'};
leg = cell(numel(folderName),1);
for i = 1:1:numel(folderName)  
    cL = crackBond(folderName{i}, 10:10:500);
    cL(cL(:,2)>0.6,:) = NaN;
    plot(cL(:,1),cL(:,2),[color(i),shape{i}],'linewidth',2)
    leg{i,1} = ['$','\it{l_b}\rm=',num2str(lcArr(i)),'m$'];
    hold on
end
xlabel('$\rm{applied\;u_2}$','interpreter','latex')
ylabel('$\rm{crack\;extend\;length(m)}$','interpreter','latex')
legend(leg,'interpreter','latex')
set(gca,'fontsize',16)